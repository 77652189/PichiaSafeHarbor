from __future__ import annotations

import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


EXPECTED_PROTECTED_AUTHORITY = {
    "slice0": {"path": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json", "size_bytes": 7175, "sha256": "f04ac9119a891885dfd7f883484f5fe1ddbfea6a1895bcd353b604dec2b294b0"},
    "slice0a": {"path": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json", "size_bytes": 2120, "sha256": "564324ad63ceb22e6728e85b09bfbcfc625ad39df3048702bc99d47c507e4b78"},
    "slice0b": {"path": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json", "size_bytes": 2864, "sha256": "83bdf169bbfa3dcfbf7b1299175fb70f27f72e3c4388672703a9df9e4c7942f3"},
}
TOOLCHAIN_FILES = (
    "src/pichia_safe_harbor/independent_transcript_source.py", "src/pichia_safe_harbor/pipeline.py",
    "scripts/run_independent_source_qualification.py", "scripts/run_independent_source_workflow.py",
    "scripts/run_test_evidence.py", "scripts/verify_independent_source.py",
    "scripts/verify_independent_source_repeatability.py", "scripts/accept_independent_source.py",
    "schemas/independent_transcript_source_manifest.v1.schema.json", "schemas/independent_transcript_source_acceptance.v1.schema.json",
    "tests/test_independent_transcript_source.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def toolchain_sha256(repo_root: Path) -> str:
    digest = hashlib.sha256()
    for relative in TOOLCHAIN_FILES:
        path = repo_root / relative
        if not path.is_file():
            raise SystemExit(f"candidate toolchain file missing: {relative}")
        digest.update(relative.encode("utf-8") + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


def verify_protected(repo_root: Path, identity: dict) -> dict:
    repo_root = repo_root.resolve()
    acceptance_path = (repo_root / identity["path"]).resolve()
    if not acceptance_path.is_relative_to(repo_root):
        raise SystemExit(f"protected acceptance escaped repository root: {identity['path']}")
    if not acceptance_path.is_file() or acceptance_path.stat().st_size != identity["size_bytes"] or sha256(acceptance_path) != identity["sha256"]:
        raise SystemExit(f"protected acceptance identity mismatch: {identity['path']}")
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    run_dir = acceptance_path.parent.resolve()
    run_manifest_path = run_dir / "run_manifest.json"
    run_identity = acceptance["run_manifest"]
    if run_manifest_path.stat().st_size != run_identity["size_bytes"] or sha256(run_manifest_path) != run_identity["sha256"]:
        raise SystemExit(f"protected run manifest mismatch: {identity['path']}")
    manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    artifact_sha256 = {}
    for name, artifact_identity in sorted(manifest["artifacts"].items()):
        path = (run_dir / name).resolve()
        if not path.is_relative_to(run_dir) or not path.is_file() or path.stat().st_size != artifact_identity["size_bytes"] or sha256(path) != artifact_identity["sha256"]:
            raise SystemExit(f"protected evidence mismatch: {identity['path']}:{name}")
        artifact_sha256[name] = artifact_identity["sha256"]
    return {"acceptance_sha256": identity["sha256"], "run_manifest_sha256": run_identity["sha256"], "artifact_sha256": artifact_sha256}


def parse_date(value: str) -> date | None:
    try:
        return date.fromisoformat(value)
    except ValueError:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    run = json.loads((args.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    reference_manifest = json.loads(args.reference_manifest.read_text(encoding="utf-8"))
    if manifest.get("protected_authority") != EXPECTED_PROTECTED_AUTHORITY:
        raise SystemExit("candidate protected authority contract mismatch")
    workflow_hash = toolchain_sha256(args.repo_root.resolve())
    if run.get("toolchain_sha256") != workflow_hash:
        raise SystemExit("candidate toolchain identity mismatch")
    if sha256(args.manifest) != run["manifest_sha256"] or sha256(args.reference_manifest) != run["reference_manifest_sha256"]:
        raise SystemExit("candidate input manifest identity mismatch")
    reference = reference_manifest["references"]["strain-b"]
    if reference["assembly_accession"] != "GCA_001746955.1":
        raise SystemExit("candidate primary reference mismatch")
    annotation_date_text = reference["annotation_accession"].rsplit(" ", 1)[-1]
    annotation_date = parse_date(annotation_date_text)
    if annotation_date is None:
        raise SystemExit("candidate annotation release date is invalid")
    qualification = json.loads((args.run_dir / "source_qualification.json").read_text(encoding="utf-8"))
    verified_artifacts = {}
    for name, identity in run["artifacts"].items():
        path = args.run_dir / name
        if path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"candidate artifact mismatch: {name}")
        verified_artifacts[name] = identity["sha256"]
    for name, identity in manifest["files"].items():
        path = args.repo_root / identity["path"]
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"candidate source identity mismatch: {name}")
    root = ET.parse(args.repo_root / manifest["files"]["sra_experiment_packages"]["path"]).getroot()
    package = next(item for item in root.findall("EXPERIMENT_PACKAGE") if item.find("./RUN_SET/RUN").get("accession") == manifest["accessions"]["run"])
    sample = package.find("SAMPLE")
    experiment = package.find("EXPERIMENT")
    run_node = package.find("./RUN_SET/RUN")
    assert sample is not None and experiment is not None and run_node is not None
    attrs = {item.findtext("TAG", ""): item.findtext("VALUE", "") for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
    normalized = next(item for item in run_node.findall("./SRAFiles/SRAFile") if item.get("semantic_name") == "SRA Normalized")
    archive = args.repo_root / manifest["files"]["sra_normalized"]["path"]
    if archive.stat().st_size != int(normalized.get("size", "0")) or md5(archive) != normalized.get("md5"):
        raise SystemExit("candidate archive declared identity mismatch")
    design = experiment.findtext("./DESIGN/DESIGN_DESCRIPTION", "")
    library = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
    assert library is not None
    library_values = [" ".join(item.itertext()).strip().lower() for item in library.iter()]
    negative_strand = {"unstranded", "not stranded", "unknown", "strand unknown", "unspecified"}
    explicit_strand = {"stranded", "directional", "forward", "reverse", "first-strand", "second-strand"}
    strand_status = "unavailable" if any(value in negative_strand for value in library_values) else "declared" if any(value in explicit_strand for value in library_values) else "unavailable"
    crossref = json.loads((args.repo_root / manifest["files"]["crossref_publication"]["path"]).read_text(encoding="utf-8"))["message"]
    license_urls = [item.get("URL", "") for item in crossref.get("license", [])]
    raw_metadata = (args.repo_root / manifest["files"]["sra_experiment_packages"]["path"]).read_text(encoding="utf-8") + (args.repo_root / manifest["files"]["ena_read_run"]["path"]).read_text(encoding="utf-8")
    wt_runs = [
        {
            "biosample": (item.find("SAMPLE").find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']").text if item.find("SAMPLE").find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']") is not None else None),
            "replicate": next((entry.findtext("VALUE", "") for entry in item.find("SAMPLE").findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE") if entry.findtext("TAG", "").lower() in {"biological_replicate", "biological replicate"}), None),
        }
        for item in root.findall("EXPERIMENT_PACKAGE")
        if item.find("SAMPLE") is not None and item.find("SAMPLE").findtext("DESCRIPTION", "") == "WT" and item.find("./RUN_SET/RUN") is not None
    ]
    biological_samples = {item["biosample"] for item in wt_runs if item["biosample"] and item["replicate"]}
    collection_date = parse_date(attrs.get("collection_date", ""))
    ena_rows = json.loads((args.repo_root / manifest["files"]["ena_read_run"]["path"]).read_text(encoding="utf-8"))
    ena_row = next((item for item in ena_rows if item.get("run_accession") == manifest["accessions"]["run"]), None)
    sample_external = sample.find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']")
    study_ref = experiment.find("STUDY_REF")
    pubmed_root = ET.parse(args.repo_root / manifest["files"]["pubmed_publication"]["path"]).getroot()
    pubmed_ids = {item.get("IdType"): (item.text or "") for item in pubmed_root.findall(".//ArticleId")}
    bioproject_root = ET.parse(args.repo_root / manifest["files"]["bioproject"]["path"]).getroot()
    bioproject_archive = bioproject_root.find(".//ArchiveID")
    crosslinks = all([
        experiment.get("accession") == manifest["accessions"]["experiment"],
        run_node.get("accession") == manifest["accessions"]["run"],
        sample_external is not None and sample_external.text == manifest["accessions"]["sample"],
        study_ref is not None and study_ref.get("accession") == manifest["accessions"]["study"],
        ena_row is not None and ena_row.get("study_accession") == manifest["accessions"]["bioproject"],
        ena_row is not None and ena_row.get("sample_accession") == manifest["accessions"]["sample"],
        ena_row is not None and ena_row.get("experiment_accession") == manifest["accessions"]["experiment"],
        bioproject_archive is not None and bioproject_archive.get("accession") == manifest["accessions"]["bioproject"],
        pubmed_ids.get("pubmed") == manifest["accessions"]["pmid"],
        pubmed_ids.get("doi") == manifest["accessions"]["publication_doi"],
    ])
    layout = next(iter(library.find("LIBRARY_LAYOUT"))).tag if library.find("LIBRARY_LAYOUT") is not None else None
    expected = {
        "archive_strain_identity": "passed" if attrs.get("strain") == "Strain-B" else "not-passed",
        "unmodified_strain-b_identity": "not-passed" if "expressing egfp" in design.lower() else "unavailable",
        "accession_crosslinks": "passed" if crosslinks else "not-passed",
        "annotation_generation_independence": "passed" if collection_date is not None and collection_date > annotation_date else "unavailable",
        "library_metadata": "passed" if {"strategy": library.findtext("LIBRARY_STRATEGY"), "source": library.findtext("LIBRARY_SOURCE"), "selection": library.findtext("LIBRARY_SELECTION"), "layout": layout} == {"strategy": "RNA-Seq", "source": "TRANSCRIPTOMIC", "selection": "cDNA", "layout": "SINGLE"} else "not-passed",
        "run_condition_specificity": "partial",
        "raw_file_license": "unresolved" if license_urls and all("springer.com/tdm" in url for url in license_urls) and "<license" not in raw_metadata.lower() and "creative commons" not in raw_metadata.lower() else "requires-review",
        "strand_specificity": strand_status,
        "biological_replication": "passed" if len(biological_samples) >= 2 else "not-passed",
        "controlled_coordinate_mapping": "not-run",
    }
    if qualification["gates"] != expected or qualification["eligible_as_independent_transcript_source"] is not False or qualification["scientific_acceptance_status"] != "blocked":
        raise SystemExit("candidate qualification gate mismatch")
    protected = {name: verify_protected(args.repo_root, identity) for name, identity in manifest["protected_authority"].items()}
    result = {
        "schema_version": 1,
        "verification_type": "independent_transcript_source_qualification",
        "status": "passed",
        "run_id": run["run_id"],
        "manifest_sha256": sha256(args.manifest),
        "reference_manifest_sha256": sha256(args.reference_manifest),
        "toolchain_sha256": workflow_hash,
        "verified_artifacts": verified_artifacts,
        "verified_source_files": {name: identity["sha256"] for name, identity in manifest["files"].items()},
        "recomputed_gates": expected,
        "eligible_as_independent_transcript_source": False,
        "scientific_acceptance_status": "blocked",
        "protected_evidence": protected,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
