from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from datetime import date, datetime
from pathlib import Path
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file, write_json
from .pipeline import _implementation_hash


EXPECTED_ARTIFACTS = {"source_qualification.json", "source_qualification.md"}
REPEATABLE_FILES = {
    *EXPECTED_ARTIFACTS,
    "run_manifest.json",
    "verification/independent_check.json",
    "verification/test_evidence.json",
}
EXPECTED_PREACCEPTANCE = {*REPEATABLE_FILES, "verification/execution_receipt.json"}
EXPECTED_PROTECTED_AUTHORITY = {
    "slice0": {"path": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json", "size_bytes": 7175, "sha256": "f04ac9119a891885dfd7f883484f5fe1ddbfea6a1895bcd353b604dec2b294b0"},
    "slice0a": {"path": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json", "size_bytes": 2120, "sha256": "564324ad63ceb22e6728e85b09bfbcfc625ad39df3048702bc99d47c507e4b78"},
    "slice0b": {"path": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json", "size_bytes": 2864, "sha256": "83bdf169bbfa3dcfbf7b1299175fb70f27f72e3c4388672703a9df9e4c7942f3"},
}
TOOLCHAIN_FILES = (
    "src/pichia_safe_harbor/independent_transcript_source.py",
    "src/pichia_safe_harbor/pipeline.py",
    "scripts/run_independent_source_qualification.py",
    "scripts/run_independent_source_workflow.py",
    "scripts/run_test_evidence.py",
    "scripts/verify_independent_source.py",
    "scripts/verify_independent_source_repeatability.py",
    "scripts/accept_independent_source.py",
    "schemas/independent_transcript_source_manifest.v1.schema.json",
    "schemas/independent_transcript_source_acceptance.v1.schema.json",
    "tests/test_independent_transcript_source.py",
)
EXPECTED_PRJNA604658_GATES = {
    "archive_strain_identity": "passed",
    "unmodified_strain-b_identity": "not-passed",
    "accession_crosslinks": "passed",
    "annotation_generation_independence": "passed",
    "library_metadata": "passed",
    "run_condition_specificity": "partial",
    "raw_file_license": "unresolved",
    "strand_specificity": "unavailable",
    "biological_replication": "not-passed",
    "controlled_coordinate_mapping": "not-run",
}


def _identity(path: Path) -> dict[str, Any]:
    return {"size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def toolchain_sha256(repo_root: Path) -> str:
    digest = hashlib.sha256()
    root = repo_root.resolve()
    for relative in TOOLCHAIN_FILES:
        path = root / relative
        if not path.is_file():
            raise ContractError(f"independent source toolchain file is missing: {relative}")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_identity(repo_root: Path, identity: dict[str, Any], label: str) -> Path:
    path = (repo_root.resolve() / identity["path"]).resolve()
    if not path.is_relative_to(repo_root.resolve()):
        raise ContractError(f"candidate source file escaped repository root: {label}")
    if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256_file(path) != identity["sha256"]:
        raise ContractError(f"candidate source file identity mismatch: {label}")
    return path


def _verify_protected_authority(repo_root: Path, identity: dict[str, Any], label: str) -> dict[str, Any]:
    acceptance_path = _verify_identity(repo_root, identity, f"protected authority {label}")
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    run_identity = acceptance.get("run_manifest")
    if not isinstance(run_identity, dict):
        raise ContractError(f"protected authority run manifest is missing: {label}")
    run_dir = acceptance_path.parent.resolve()
    run_manifest_path = run_dir / "run_manifest.json"
    if not run_manifest_path.is_file() or _identity(run_manifest_path) != run_identity:
        raise ContractError(f"protected authority run manifest mismatch: {label}")
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    artifacts = run_manifest.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ContractError(f"protected authority artifacts are missing: {label}")
    artifact_sha256 = {}
    for name, artifact_identity in sorted(artifacts.items()):
        artifact_path = (run_dir / name).resolve()
        if not artifact_path.is_relative_to(run_dir) or not artifact_path.is_file() or _identity(artifact_path) != artifact_identity:
            raise ContractError(f"protected authority artifact mismatch: {label}:{name}")
        artifact_sha256[name] = artifact_identity["sha256"]
    return {
        "acceptance_sha256": identity["sha256"],
        "run_manifest_sha256": run_identity["sha256"],
        "artifact_sha256": artifact_sha256,
    }


def _selected_package(root: ET.Element, run_accession: str) -> ET.Element:
    for package in root.findall("EXPERIMENT_PACKAGE"):
        run = package.find("./RUN_SET/RUN")
        if run is not None and run.get("accession") == run_accession:
            return package
    raise ContractError(f"candidate run missing from SRA XML: {run_accession}")


def _iso_date(value: str) -> date | None:
    try:
        return date.fromisoformat(value)
    except ValueError:
        return None


def _strand_status(library: ET.Element) -> str:
    values = [" ".join(item.itertext()).strip().lower() for item in library.iter()]
    negatives = {"unstranded", "not stranded", "unknown", "strand unknown", "unspecified"}
    if any(value in negatives for value in values):
        return "unavailable"
    explicit = {"stranded", "directional", "forward", "reverse", "first-strand", "second-strand"}
    return "declared" if any(value in explicit for value in values) else "unavailable"


def qualify_independent_transcript_source(
    manifest_path: Path,
    reference_manifest_path: Path,
    output_dir: Path,
    repo_root: Path,
) -> dict[str, Any]:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise ContractError(f"output directory already exists: {output_dir}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("candidate_id") != "prjna604658-srr11011828":
        raise ContractError("candidate manifest identity mismatch")
    if manifest.get("protected_authority") != EXPECTED_PROTECTED_AUTHORITY:
        raise ContractError("candidate protected authority contract mismatch")
    files = {name: _verify_identity(repo_root, identity, name) for name, identity in manifest["files"].items()}
    sra_root = ET.parse(files["sra_experiment_packages"]).getroot()
    package = _selected_package(sra_root, manifest["accessions"]["run"])
    experiment = package.find("EXPERIMENT")
    sample = package.find("SAMPLE")
    run = package.find("./RUN_SET/RUN")
    assert experiment is not None and sample is not None and run is not None
    ena_rows = json.loads(files["ena_read_run"].read_text(encoding="utf-8"))
    ena_row = next((item for item in ena_rows if item.get("run_accession") == manifest["accessions"]["run"]), None)
    pubmed_root = ET.parse(files["pubmed_publication"]).getroot()
    bioproject_root = ET.parse(files["bioproject"]).getroot()
    sample_external = sample.find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']")
    study_ref = experiment.find("STUDY_REF")
    pubmed_ids = {item.get("IdType"): (item.text or "") for item in pubmed_root.findall(".//ArticleId")}
    bioproject_archive = bioproject_root.find(".//ArchiveID")
    accession_crosslinks = all([
        experiment.get("accession") == manifest["accessions"]["experiment"],
        run.get("accession") == manifest["accessions"]["run"],
        sample_external is not None and sample_external.text == manifest["accessions"]["sample"],
        study_ref is not None and study_ref.get("accession") == manifest["accessions"]["study"],
        ena_row is not None and ena_row.get("study_accession") == manifest["accessions"]["bioproject"],
        ena_row is not None and ena_row.get("experiment_accession") == manifest["accessions"]["experiment"],
        ena_row is not None and ena_row.get("sample_accession") == manifest["accessions"]["sample"],
        bioproject_archive is not None and bioproject_archive.get("accession") == manifest["accessions"]["bioproject"],
        pubmed_ids.get("pubmed") == manifest["accessions"]["pmid"],
        pubmed_ids.get("doi") == manifest["accessions"]["publication_doi"],
    ])
    if not accession_crosslinks:
        raise ContractError("candidate accession cross-file identity mismatch")
    attributes = {item.findtext("TAG", ""): item.findtext("VALUE", "") for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
    library = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
    assert library is not None
    sra_file = next((item for item in run.findall("./SRAFiles/SRAFile") if item.get("semantic_name") == "SRA Normalized"), None)
    if sra_file is None:
        raise ContractError("candidate normalized SRA identity is missing")
    archive_path = files["sra_normalized"]
    declared_md5 = sra_file.get("md5")
    actual_md5 = _md5_file(archive_path)
    if archive_path.stat().st_size != int(sra_file.get("size", "0")) or actual_md5 != declared_md5:
        raise ContractError("candidate normalized SRA size or MD5 mismatch")
    crossref = json.loads(files["crossref_publication"].read_text(encoding="utf-8"))["message"]
    if crossref.get("DOI") != manifest["accessions"]["publication_doi"]:
        raise ContractError("candidate Crossref DOI mismatch")
    license_urls = sorted(item.get("URL", "") for item in crossref.get("license", []))
    raw_metadata = files["sra_experiment_packages"].read_text(encoding="utf-8") + files["ena_read_run"].read_text(encoding="utf-8")
    raw_file_license = "unresolved" if "creative commons" not in raw_metadata.lower() and "<license" not in raw_metadata.lower() else "requires-review"
    reference_manifest = json.loads(reference_manifest_path.read_text(encoding="utf-8"))
    reference_record = reference_manifest["references"]["strain-b"]
    if reference_record.get("assembly_accession") != "GCA_001746955.1":
        raise ContractError("candidate qualification changed the locked primary assembly")
    annotation_accession = reference_record["annotation_accession"]
    annotation_date = annotation_accession.rsplit(" ", 1)[-1]
    collection_date = attributes.get("collection_date", "")
    collection_value = _iso_date(collection_date)
    annotation_value = _iso_date(annotation_date)
    independence = "passed" if collection_value is not None and annotation_value is not None and collection_value > annotation_value else "unavailable"
    design = experiment.findtext("./DESIGN/DESIGN_DESCRIPTION", "")
    strand_status = _strand_status(library)
    wt_records = [
        {
            "run": p.find("./RUN_SET/RUN").get("accession"),
            "biosample": (p.find("SAMPLE").find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']").text if p.find("SAMPLE").find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']") is not None else None),
            "replicate": next((item.findtext("VALUE", "") for item in p.find("SAMPLE").findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE") if item.findtext("TAG", "").lower() in {"biological_replicate", "biological replicate"}), None),
        }
        for p in sra_root.findall("EXPERIMENT_PACKAGE")
        if p.find("SAMPLE") is not None and p.find("SAMPLE").findtext("DESCRIPTION", "") == "WT" and p.find("./RUN_SET/RUN") is not None
    ]
    biological_samples = {item["biosample"] for item in wt_records if item["biosample"] and item["replicate"]}
    replicate_status = "passed" if len(biological_samples) >= 2 else "not-passed"
    library_values = {
        "strategy": library.findtext("LIBRARY_STRATEGY"),
        "source": library.findtext("LIBRARY_SOURCE"),
        "selection": library.findtext("LIBRARY_SELECTION"),
        "layout": next(iter(library.find("LIBRARY_LAYOUT"))).tag if library.find("LIBRARY_LAYOUT") is not None else None,
    }
    library_gate = "passed" if library_values == {"strategy": "RNA-Seq", "source": "TRANSCRIPTOMIC", "selection": "cDNA", "layout": "SINGLE"} else "not-passed"
    gates = {
        "archive_strain_identity": "passed" if attributes.get("strain") == "Strain-B" else "not-passed",
        "unmodified_strain-b_identity": "not-passed" if "expressing egfp" in design.lower() else "unavailable",
        "accession_crosslinks": "passed",
        "annotation_generation_independence": independence,
        "library_metadata": library_gate,
        "run_condition_specificity": "partial",
        "raw_file_license": raw_file_license,
        "strand_specificity": strand_status,
        "biological_replication": replicate_status,
        "controlled_coordinate_mapping": "not-run",
    }
    if gates != EXPECTED_PRJNA604658_GATES:
        raise ContractError("candidate qualification gates differ from the reviewed PRJNA604658 contract")
    eligible = all(value == "passed" for value in gates.values())
    qualification = {
        "schema_version": 1,
        "candidate_id": manifest["candidate_id"],
        "accessions": manifest["accessions"],
        "target_coordinate_space": "GCA_001746955.1",
        "exact_target_strain_coordinates": False,
        "archive_identity": {
            "strain": attributes.get("strain"),
            "sample_description": sample.findtext("DESCRIPTION"),
            "design_description": design,
            "collection_date": collection_date,
            "sample_type": attributes.get("sample_type"),
        },
        "experiment": {
            **library_values,
            "instrument": experiment.findtext("./PLATFORM/ILLUMINA/INSTRUMENT_MODEL"),
            "read_count": int(run.get("total_spots", "0")),
            "base_count": int(run.get("total_bases", "0")),
            "wt_records": wt_records,
            "study_context": "methanol induction for 120 hours is stated in the study abstract but is not bound to the selected run metadata",
            "selected_run_condition": design,
        },
        "file_identity": {
            "path": manifest["files"]["sra_normalized"]["path"],
            "size_bytes": archive_path.stat().st_size,
            "md5": actual_md5,
            "sha256": sha256_file(archive_path),
            "ncbi_declared_md5": declared_md5,
        },
        "publication": {
            "doi": crossref.get("DOI"),
            "published": crossref.get("published", {}).get("date-parts", [[None]])[0],
            "license_urls": license_urls,
            "license_interpretation": "Springer TDM terms apply to article use and are not promoted to a raw sequencing file license",
        },
        "generation_independence_evidence": {
            "current_annotation_release": annotation_date,
            "sample_collection_date": collection_date,
            "interpretation": "the 2019 sample cannot have participated in the 2016 assembly/annotation generation process",
        },
        "gates": gates,
        "eligible_as_independent_transcript_source": eligible,
        "scientific_acceptance_status": "accepted" if eligible else "blocked",
        "stop_line": "do not start full-genome reannotation or Slice 1",
    }
    implementation_hash = _implementation_hash()
    workflow_hash = toolchain_sha256(Path(__file__).resolve().parents[2])
    material = json.dumps({
        "manifest_sha256": sha256_file(manifest_path),
        "reference_manifest_sha256": sha256_file(reference_manifest_path),
        "implementation_sha256": implementation_hash,
        "toolchain_sha256": workflow_hash,
    }, sort_keys=True).encode("utf-8")
    run_id = "independent-source-" + hashlib.sha256(material).hexdigest()[:16]
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        write_json(temp_dir / "source_qualification.json", qualification)
        report = "\n".join([
            "# Independent transcript source qualification: PRJNA604658 / SRR11011828",
            "",
            "## Result",
            "",
            "- Engineering qualification completed; scientific status remains blocked.",
            "- The archive labels the selected sample Strain-B, but the experiment is WT expressing eGFP and must not be described as unmodified Strain-B.",
            "- The 2019 sample is temporally independent of the 2016 GCA_001746955.1 assembly/annotation generation process.",
            "- The raw sequencing file has no explicit per-file license; Springer TDM terms are not a raw-data license.",
            "- Library strand specificity is unavailable and the WT condition has one run only.",
            "- No controlled mapping or four-chromosome quality probe was run because metadata qualification already failed.",
            "",
            "## Stop line",
            "",
            "Do not start full-genome reannotation, candidate generation, threshold freezing, or Slice 1.",
        ]) + "\n"
        (temp_dir / "source_qualification.md").write_text(report, encoding="utf-8")
        artifacts = {name: _identity(temp_dir / name) for name in sorted(EXPECTED_ARTIFACTS)}
        run_manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "execution_status": "complete",
            "verification_status": "not_run",
            "scientific_acceptance_status": "blocked",
            "candidate_id": manifest["candidate_id"],
            "manifest_sha256": sha256_file(manifest_path),
            "reference_manifest_sha256": sha256_file(reference_manifest_path),
            "implementation_sha256": implementation_hash,
            "toolchain_sha256": workflow_hash,
            "source_files": {name: _identity(path) | {"path": manifest["files"][name]["path"]} for name, path in files.items()},
            "artifacts": artifacts,
        }
        write_json(temp_dir / "run_manifest.json", run_manifest)
        os.replace(temp_dir, output_dir)
        return run_manifest
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


def create_independent_source_acceptance(
    run_dir: Path,
    peer_run_dir: Path,
    independent_path: Path,
    test_evidence_path: Path,
    repeatability_path: Path,
    repo_root: Path,
    manifest_path: Path,
    reference_manifest_path: Path,
) -> dict[str, Any]:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    qualification = json.loads((run_dir / "source_qualification.json").read_text(encoding="utf-8"))
    independent = json.loads(independent_path.read_text(encoding="utf-8"))
    tests = json.loads(test_evidence_path.read_text(encoding="utf-8"))
    repeatability = json.loads(repeatability_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    workflow_hash = toolchain_sha256(Path(__file__).resolve().parents[2])
    if manifest.get("protected_authority") != EXPECTED_PROTECTED_AUTHORITY:
        raise ContractError("independent source protected authority contract mismatch")
    if run.get("implementation_sha256") != _implementation_hash() or run.get("toolchain_sha256") != workflow_hash or run.get("manifest_sha256") != sha256_file(manifest_path) or run.get("reference_manifest_sha256") != sha256_file(reference_manifest_path):
        raise ContractError("independent source acceptance input identity mismatch")
    expected_source_files = {}
    for name, identity in manifest["files"].items():
        path = _verify_identity(repo_root, identity, name)
        expected_source_files[name] = _identity(path) | {"path": identity["path"]}
    if run.get("source_files") != expected_source_files:
        raise ContractError("independent source acceptance source identity mismatch")
    if set(run.get("artifacts", {})) != EXPECTED_ARTIFACTS or qualification.get("eligible_as_independent_transcript_source") is not False or qualification.get("scientific_acceptance_status") != "blocked":
        raise ContractError("independent source eligibility was overclaimed")
    if qualification.get("gates") != EXPECTED_PRJNA604658_GATES:
        raise ContractError("independent source gate summary mismatch")
    for name, identity in run["artifacts"].items():
        path = run_dir / name
        if not path.is_file() or _identity(path) != identity:
            raise ContractError(f"independent source artifact mismatch: {name}")
    expected_hashes = {name: identity["sha256"] for name, identity in run["artifacts"].items()}
    expected_source_hashes = {name: identity["sha256"] for name, identity in run["source_files"].items()}
    if independent_path.resolve() != (run_dir / "verification/independent_check.json").resolve() or test_evidence_path.resolve() != (run_dir / "verification/test_evidence.json").resolve():
        raise ContractError("independent source verification paths are not bound to the accepted run")
    if independent.get("verification_type") != "independent_transcript_source_qualification" or independent.get("status") != "passed" or independent.get("run_id") != run["run_id"] or independent.get("manifest_sha256") != run["manifest_sha256"] or independent.get("reference_manifest_sha256") != run["reference_manifest_sha256"] or independent.get("toolchain_sha256") != workflow_hash or independent.get("verified_artifacts") != expected_hashes or independent.get("verified_source_files") != expected_source_hashes or independent.get("recomputed_gates") != EXPECTED_PRJNA604658_GATES or independent.get("eligible_as_independent_transcript_source") is not False or independent.get("scientific_acceptance_status") != "blocked":
        raise ContractError("independent source verification mismatch")
    protected = independent.get("protected_evidence", {})
    expected_protected = {
        name: _verify_protected_authority(repo_root, identity, name)
        for name, identity in manifest["protected_authority"].items()
    }
    if protected != expected_protected:
        raise ContractError("independent source protected evidence mismatch")
    if tests.get("schema_version") != 1 or tests.get("evidence_type") != "automated_tests" or tests.get("status") != "passed" or tests.get("command") != ["python", "-m", "pytest", "-q"] or tests.get("passed_count", 0) < 70 or tests.get("implementation_sha256") != run["implementation_sha256"] or tests.get("toolchain_sha256") != workflow_hash:
        raise ContractError("independent source test evidence mismatch")
    if repeatability.get("verification_type") != "independent_transcript_source_repeatability" or repeatability.get("status") != "passed" or repeatability.get("run_id") != run["run_id"] or repeatability.get("toolchain_sha256") != workflow_hash or repeatability.get("all_files_identical") is not True or set(repeatability.get("file_sha256", {})) != REPEATABLE_FILES:
        raise ContractError("independent source repeatability mismatch")
    if run_dir.resolve() == peer_run_dir.resolve() or set(repeatability.get("run_directories", [])) != {run_dir.name, peer_run_dir.name}:
        raise ContractError("independent source peer identity mismatch")
    run_files = {path.relative_to(run_dir).as_posix() for path in run_dir.rglob("*") if path.is_file()}
    peer_files = {path.relative_to(peer_run_dir).as_posix() for path in peer_run_dir.rglob("*") if path.is_file()}
    if run_files != EXPECTED_PREACCEPTANCE or peer_files != EXPECTED_PREACCEPTANCE:
        raise ContractError("independent source run file set mismatch")
    for relative, digest in repeatability["file_sha256"].items():
        if sha256_file(run_dir / relative) != digest or sha256_file(peer_run_dir / relative) != digest:
            raise ContractError(f"independent source repeatability hash mismatch: {relative}")
    if repeatability["file_sha256"]["verification/independent_check.json"] != sha256_file(independent_path) or repeatability["file_sha256"]["verification/test_evidence.json"] != sha256_file(test_evidence_path):
        raise ContractError("independent source verification evidence is not the repeated evidence")
    receipts = repeatability.get("execution_receipts", {})
    if set(receipts) != {run_dir.name, peer_run_dir.name}:
        raise ContractError("independent source execution receipt set mismatch")
    invocation_ids = set()
    coordinator_ids = set()
    run_slots = set()
    for directory, base in ((run_dir.name, run_dir), (peer_run_dir.name, peer_run_dir)):
        receipt_path = base / "verification/execution_receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        expected_receipt = receipts[directory]
        steps = receipt.get("workflow_steps", [])
        code_root = Path(__file__).resolve().parents[2]
        expected_commands = [
            [sys.executable, str(code_root / "scripts/run_independent_source_qualification.py"), "--manifest", str(manifest_path.resolve()), "--reference-manifest", str(reference_manifest_path.resolve()), "--output-dir", str(base.resolve()), "--repo-root", str(repo_root.resolve())],
            [sys.executable, str(code_root / "scripts/run_test_evidence.py"), "--output", str((base / "verification/test_evidence.json").resolve())],
            [sys.executable, str(code_root / "scripts/verify_independent_source.py"), "--run-dir", str(base.resolve()), "--manifest", str(manifest_path.resolve()), "--reference-manifest", str(reference_manifest_path.resolve()), "--repo-root", str(repo_root.resolve()), "--output", str((base / "verification/independent_check.json").resolve())],
        ]
        commands_valid = len(steps) == 3 and all(step.get("command") == expected for step, expected in zip(steps, expected_commands))
        times_valid = True
        try:
            intervals = [(datetime.fromisoformat(step["started_at_utc"]), datetime.fromisoformat(step["completed_at_utc"])) for step in steps]
            times_valid = all(start <= end for start, end in intervals) and all(intervals[index][1] <= intervals[index + 1][0] for index in range(len(intervals) - 1))
        except (KeyError, TypeError, ValueError):
            times_valid = False
        outputs_valid = receipt.get("workflow_outputs") == {
            relative: _identity(base / relative) for relative in sorted(REPEATABLE_FILES)
        }
        if (
            _identity(receipt_path) != expected_receipt["identity"]
            or receipt.get("evidence_type") != "independent_execution_receipt"
            or receipt.get("run_id") != run["run_id"]
            or receipt.get("candidate_id") != run["candidate_id"]
            or receipt.get("invocation_id") != expected_receipt["invocation_id"]
            or receipt.get("manifest_sha256") != run["manifest_sha256"]
            or receipt.get("reference_manifest_sha256") != run["reference_manifest_sha256"]
            or receipt.get("implementation_sha256") != run["implementation_sha256"]
            or receipt.get("toolchain_sha256") != workflow_hash
            or not isinstance(receipt.get("command"), list)
            or "run_independent_source_qualification.py" not in " ".join(receipt["command"])
            or receipt.get("workflow_evidence_type") != "orchestrated_independent_source_execution"
            or receipt.get("workflow_status") != "passed"
            or receipt.get("workflow_cwd") != str(code_root)
            or not commands_valid
            or receipt["command"][1:] != steps[0]["command"][1:]
            or any(step.get("exit_code") != 0 or not step.get("stdout_sha256") for step in steps)
            or not times_valid
            or not outputs_valid
        ):
            raise ContractError("independent source execution receipt mismatch")
        invocation_ids.add(receipt["invocation_id"])
        coordinator_ids.add(receipt.get("coordinator_id"))
        run_slots.add(receipt.get("run_slot"))
    if len(invocation_ids) != 2 or len(coordinator_ids) != 1 or None in coordinator_ids or run_slots != {"run-a", "run-b"}:
        raise ContractError("independent source repeat runs do not have distinct execution receipts")
    result = {
        "schema_version": 1,
        "acceptance_version": "independent-transcript-source-v1",
        "run_id": run["run_id"],
        "candidate_id": run["candidate_id"],
        "toolchain_sha256": workflow_hash,
        "execution_status": "complete",
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
        "eligible_as_independent_transcript_source": False,
        "qualification_gates": qualification["gates"],
        "blocking_gates": {key: value for key, value in EXPECTED_PRJNA604658_GATES.items() if value != "passed"},
        "run_manifest": _identity(run_dir / "run_manifest.json"),
        "artifacts": run["artifacts"],
        "verification_evidence": {
            "independent": _identity(independent_path),
            "automated_tests": _identity(test_evidence_path) | {"passed_count": tests["passed_count"]},
            "repeatability": _identity(repeatability_path),
        },
        "stop_line": "remain blocked; do not start full-genome reannotation or Slice 1",
    }
    write_json(run_dir / "acceptance_manifest.json", result)
    return result
