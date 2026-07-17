from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
import uuid
import xml.etree.ElementTree as ET
from datetime import date, datetime, timezone
from pathlib import Path


CANDIDATE_ID = "prjna1210090-srp557139-wt-no-stress"
SELECTED = {
    "SRR31989016": {"experiment": "SRX27343955", "sample": "SAMN46238833", "replicate": "3"},
    "SRR31989027": {"experiment": "SRX27343944", "sample": "SAMN46238832", "replicate": "2"},
    "SRR31989028": {"experiment": "SRX27343943", "sample": "SAMN46238831", "replicate": "1"},
}
GATES = {
    "archive_strain_identity": "passed",
    "unmodified_strain-b_identity": "unresolved",
    "accession_crosslinks": "passed",
    "annotation_generation_independence": "passed",
    "library_metadata": "passed",
    "run_condition_specificity": "passed",
    "raw_file_license": "unresolved",
    "strand_specificity": "unavailable",
    "biological_replication": "passed",
    "raw_file_identity": "not-acquired",
    "controlled_coordinate_mapping": "not-run",
}
EXPECTED_PROTECTED = {
    "prjna604658": {"path": "local_runs/independent_transcript_sources/prjna604658/qualification_v6_run1/acceptance_manifest.json", "sha256": "d2af88922a46945a215fe8eaf7bf9b724d68afb0dd539b283e5fe02ba9db0a47", "size_bytes": 2255},
    "slice0": {"path": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json", "sha256": "f04ac9119a891885dfd7f883484f5fe1ddbfea6a1895bcd353b604dec2b294b0", "size_bytes": 7175},
    "slice0a": {"path": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json", "sha256": "564324ad63ceb22e6728e85b09bfbcfc625ad39df3048702bc99d47c507e4b78", "size_bytes": 2120},
    "slice0b": {"path": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json", "sha256": "83bdf169bbfa3dcfbf7b1299175fb70f27f72e3c4388672703a9df9e4c7942f3", "size_bytes": 2864},
}
TOOLCHAIN = (
    "scripts/acquire_prjna1210090_metadata.py",
    "scripts/run_prjna1210090_metadata.py",
    "scripts/verify_prjna1210090_metadata.py",
    "scripts/run_prjna1210090_workflow.py",
    "schemas/prjna1210090_metadata_manifest.v1.schema.json",
    "schemas/prjna1210090_metadata_acceptance.v1.schema.json",
    "tests/test_prjna1210090_metadata.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": sha256(path)}


def toolchain_sha256(repo_root: Path) -> str:
    digest = hashlib.sha256()
    for relative in TOOLCHAIN:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"toolchain file missing: {relative}")
        digest.update(relative.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


def verified_file(repo_root: Path, value: dict, label: str) -> Path:
    root = repo_root.resolve()
    path = (root / value["path"]).resolve()
    if not path.is_relative_to(root) or not path.is_file() or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}:
        raise ValueError(f"file identity mismatch: {label}")
    return path


def verify_protected(repo_root: Path, value: dict, label: str) -> dict:
    acceptance_path = verified_file(repo_root, value, label)
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    run_identity = acceptance.get("run_manifest")
    if not isinstance(run_identity, dict):
        raise ValueError(f"protected run manifest missing: {label}")
    run_manifest_path = acceptance_path.parent / "run_manifest.json"
    if identity(run_manifest_path) != run_identity:
        raise ValueError(f"protected run manifest mismatch: {label}")
    run = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    artifact_sha256 = {}
    for name, artifact_identity in sorted(run["artifacts"].items()):
        path = (acceptance_path.parent / name).resolve()
        if not path.is_relative_to(acceptance_path.parent.resolve()) or identity(path) != artifact_identity:
            raise ValueError(f"protected artifact mismatch: {label}:{name}")
        artifact_sha256[name] = artifact_identity["sha256"]
    return {"acceptance_sha256": value["sha256"], "run_manifest_sha256": run_identity["sha256"], "artifact_sha256": artifact_sha256}


def qualify(manifest_path: Path, output_dir: Path, repo_root: Path) -> dict:
    if output_dir.exists():
        raise ValueError(f"output exists: {output_dir}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("candidate_id") != CANDIDATE_ID or manifest.get("accessions", {}).get("runs") != SELECTED:
        raise ValueError("candidate manifest contract mismatch")
    if manifest.get("protected_authority") != EXPECTED_PROTECTED:
        raise ValueError("protected authority contract mismatch")
    protected_evidence = {name: verify_protected(repo_root, value, name) for name, value in manifest["protected_authority"].items()}
    files = {name: verified_file(repo_root, value, name) for name, value in manifest["files"].items()}
    reference = json.loads(files["reference_manifest"].read_text(encoding="utf-8"))["references"]["strain-b"]
    if reference.get("assembly_accession") != "GCA_001746955.1":
        raise ValueError("locked reference assembly mismatch")
    reference_date = reference["annotation_accession"].rsplit(" ", 1)[-1]
    if reference_date != manifest["reference_annotation_release"]:
        raise ValueError("reference annotation release mismatch")
    acquisition = json.loads(files["acquisition_evidence"].read_text(encoding="utf-8"))
    if acquisition.get("raw_read_acquisition", {}).get("status") != "not-performed" or acquisition.get("selected_runs") != SELECTED:
        raise ValueError("acquisition evidence contract mismatch")
    sra = ET.parse(files["sra_experiment_packages"]).getroot()
    ena_rows = {row["run_accession"]: row for row in json.loads(files["ena_read_run"].read_text(encoding="utf-8"))}
    bioproject = ET.parse(files["bioproject"]).getroot()
    archive_id = bioproject.find(".//ArchiveID")
    bioproject_text = " ".join(bioproject.itertext()).lower()
    if archive_id is None or archive_id.get("accession") != manifest["accessions"]["bioproject"]:
        raise ValueError("BioProject accession mismatch")
    if "engineered k. phaffii" not in bioproject_text or "cell factory" not in bioproject_text:
        raise ValueError("BioProject engineering context mismatch")
    records = []
    for run_accession, expected in SELECTED.items():
        package = next((item for item in sra.findall(".//EXPERIMENT_PACKAGE") if item.find("./RUN_SET/RUN") is not None and item.find("./RUN_SET/RUN").get("accession") == run_accession), None)
        if package is None:
            raise ValueError(f"selected run missing: {run_accession}")
        experiment = package.find("EXPERIMENT")
        sample = package.find("SAMPLE")
        library = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
        attrs = {item.findtext("TAG", ""): item.findtext("VALUE", "") for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
        biosample = next((item.text for item in sample.findall("./IDENTIFIERS/EXTERNAL_ID") if item.get("namespace") == "BioSample"), None)
        layout = next(iter(library.find("LIBRARY_LAYOUT"))).tag
        design = experiment.findtext("./DESIGN/DESIGN_DESCRIPTION", "")
        ena = ena_rows[run_accession]
        record = {
            "run": run_accession, "experiment": experiment.get("accession"), "sample": biosample,
            "replicate": attrs.get("replicate", "").removeprefix("biological replicate "),
            "strain": attrs.get("strain"), "collection_date": attrs.get("collection_date"), "design": design,
            "library": {"strategy": library.findtext("LIBRARY_STRATEGY"), "source": library.findtext("LIBRARY_SOURCE"), "selection": library.findtext("LIBRARY_SELECTION"), "layout": layout, "instrument": ena.get("instrument_model")},
            "declared_fastq": {"bytes": ena.get("fastq_bytes", "").split(";"), "md5": ena.get("fastq_md5", "").split(";"), "ftp": ena.get("fastq_ftp", "").split(";")},
        }
        if record["experiment"] != expected["experiment"] or record["sample"] != expected["sample"] or record["replicate"] != expected["replicate"] or ena.get("experiment_accession") != expected["experiment"] or ena.get("sample_accession") != expected["sample"] or ena.get("study_accession") != manifest["accessions"]["bioproject"] or ena.get("secondary_study_accession") != manifest["accessions"]["study"]:
            raise ValueError(f"accession crosslink mismatch: {run_accession}")
        if record["strain"] != "GS115" or "wild type" not in design.lower() or "without ferulic acid stress" not in design.lower():
            raise ValueError(f"selected run identity mismatch: {run_accession}")
        if record["library"] != {"strategy": "RNA-Seq", "source": "TRANSCRIPTOMIC", "selection": "RANDOM PCR", "layout": "PAIRED", "instrument": "BGISEQ-500"}:
            raise ValueError(f"library mismatch: {run_accession}")
        if date.fromisoformat(record["collection_date"]) <= date.fromisoformat(reference_date):
            raise ValueError(f"annotation independence mismatch: {run_accession}")
        records.append(record)
    if {item["replicate"] for item in records} != {"1", "2", "3"} or len({item["sample"] for item in records}) != 3:
        raise ValueError("biological replication mismatch")
    metadata = (files["sra_experiment_packages"].read_text(encoding="utf-8") + files["ena_read_run"].read_text(encoding="utf-8")).lower()
    policy = (files["ncbi_policy"].read_text(encoding="utf-8", errors="ignore") + files["ena_terms"].read_text(encoding="utf-8", errors="ignore")).lower()
    if "no rights to transfer to a third party" not in policy or "original data may be subject to rights claimed by third parties" not in policy:
        raise ValueError("license policy evidence mismatch")
    if any(token in metadata for token in ("first-strand", "second-strand", "strand-specific", "stranded")):
        raise ValueError("unexpected strand declaration")
    workflow_hash = toolchain_sha256(Path(__file__).resolve().parents[1])
    material = json.dumps({"manifest_sha256": sha256(manifest_path), "toolchain_sha256": workflow_hash}, sort_keys=True).encode()
    run_id = "metadata-source-" + hashlib.sha256(material).hexdigest()[:16]
    result = {
        "schema_version": 1, "candidate_id": CANDIDATE_ID, "accessions": manifest["accessions"],
        "selected_condition": "wild-type GS115 without ferulic acid stress",
        "strain_identity_interpretation": "selected runs are labelled wild type GS115, but the BioProject describes an engineered vanillin-producing cell factory; unmodified GS115 identity remains unresolved",
        "records": records,
        "license_evidence": {"status": "unresolved", "interpretation": "NCBI/ENA access and repository policies do not grant a run-specific or FASTQ-specific license"},
        "strand_specificity": {"status": "unavailable", "interpretation": "RANDOM PCR and paired layout do not establish first/second-strand orientation"},
        "raw_reads": {"status": "not-acquired", "declared_total_bytes": sum(int(value) for item in records for value in item["declared_fastq"]["bytes"])},
        "gates": GATES, "eligible_as_independent_transcript_source": False,
        "scientific_acceptance_status": "blocked",
        "stop_line": "metadata hard gates failed; do not download raw reads, map coordinates, reannotate, or enter Slice 1",
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        (temp / "metadata_qualification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (temp / "metadata_qualification.md").write_text("# PRJNA1210090 metadata qualification\n\n- Selected runs are labelled wild-type GS115, but the BioProject describes an engineered vanillin-producing cell factory; unmodified GS115 identity remains unresolved.\n- Condition binding, temporal independence, library metadata, and three biological replicates passed.\n- Strand orientation is unavailable and raw FASTQ licensing is unresolved.\n- Raw reads were not downloaded; controlled mapping was not run.\n- Scientific status remains blocked; do not enter Slice 1.\n", encoding="utf-8")
        artifacts = {name: identity(temp / name) for name in ("metadata_qualification.json", "metadata_qualification.md")}
        run_manifest = {
            "schema_version": 1, "run_id": run_id, "candidate_id": CANDIDATE_ID,
            "execution_status": "complete", "verification_status": "not_run", "scientific_acceptance_status": "blocked",
            "manifest_sha256": sha256(manifest_path), "toolchain_sha256": workflow_hash,
            "protected_evidence": protected_evidence,
            "source_files": {name: identity(path) | {"path": manifest["files"][name]["path"]} for name, path in files.items()},
            "artifacts": artifacts,
        }
        (temp / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temp, output_dir)
        receipt = {
            "schema_version": 1, "evidence_type": "independent_execution_receipt", "candidate_id": CANDIDATE_ID,
            "run_id": run_id, "invocation_id": str(uuid.uuid4()), "producer_completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "manifest_sha256": run_manifest["manifest_sha256"], "toolchain_sha256": workflow_hash,
        }
        receipt_path = output_dir / "verification/execution_receipt.json"
        receipt_path.parent.mkdir(parents=True)
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return run_manifest
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    print(json.dumps(qualify(args.manifest, args.output_dir, args.repo_root), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
