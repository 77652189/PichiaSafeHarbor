from __future__ import annotations

import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


CANDIDATE_ID = "prjna1210090-srp557139-wt-no-stress"
SELECTED = {
    "SRR31989016": {"experiment": "SRX27343955", "sample": "SAMN46238833", "replicate": "3"},
    "SRR31989027": {"experiment": "SRX27343944", "sample": "SAMN46238832", "replicate": "2"},
    "SRR31989028": {"experiment": "SRX27343943", "sample": "SAMN46238831", "replicate": "1"},
}
EXPECTED_GATES = {
    "archive_strain_identity": "passed", "unmodified_strain-b_identity": "unresolved", "accession_crosslinks": "passed",
    "annotation_generation_independence": "passed", "library_metadata": "passed", "run_condition_specificity": "passed",
    "raw_file_license": "unresolved", "strand_specificity": "unavailable", "biological_replication": "passed",
    "raw_file_identity": "not-acquired", "controlled_coordinate_mapping": "not-run",
}
EXPECTED_PROTECTED = {
    "prjna604658": {"path": "local_runs/independent_transcript_sources/prjna604658/qualification_v6_run1/acceptance_manifest.json", "sha256": "d2af88922a46945a215fe8eaf7bf9b724d68afb0dd539b283e5fe02ba9db0a47", "size_bytes": 2255},
    "slice0": {"path": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json", "sha256": "f04ac9119a891885dfd7f883484f5fe1ddbfea6a1895bcd353b604dec2b294b0", "size_bytes": 7175},
    "slice0a": {"path": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json", "sha256": "564324ad63ceb22e6728e85b09bfbcfc625ad39df3048702bc99d47c507e4b78", "size_bytes": 2120},
    "slice0b": {"path": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json", "sha256": "83bdf169bbfa3dcfbf7b1299175fb70f27f72e3c4388672703a9df9e4c7942f3", "size_bytes": 2864},
}
TOOLCHAIN = (
    "scripts/acquire_prjna1210090_metadata.py", "scripts/run_prjna1210090_metadata.py",
    "scripts/verify_prjna1210090_metadata.py", "scripts/run_prjna1210090_workflow.py",
    "schemas/prjna1210090_metadata_manifest.v1.schema.json", "schemas/prjna1210090_metadata_acceptance.v1.schema.json",
    "tests/test_prjna1210090_metadata.py",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": digest(path)}


def toolchain(root: Path) -> str:
    value = hashlib.sha256()
    for relative in TOOLCHAIN:
        path = root / relative
        value.update(relative.encode() + b"\0" + path.read_bytes() + b"\0")
    return value.hexdigest()


def verify_protected(root: Path, value: dict, label: str) -> dict:
    acceptance_path = (root / value["path"]).resolve()
    if not acceptance_path.is_relative_to(root) or identity(acceptance_path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}:
        raise SystemExit(f"protected acceptance mismatch: {label}")
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8")); run_identity = acceptance["run_manifest"]
    run_path = acceptance_path.parent / "run_manifest.json"
    if identity(run_path) != run_identity:
        raise SystemExit(f"protected run mismatch: {label}")
    protected_run = json.loads(run_path.read_text(encoding="utf-8")); artifacts = {}
    for name, artifact_identity in sorted(protected_run["artifacts"].items()):
        path = (acceptance_path.parent / name).resolve()
        if not path.is_relative_to(acceptance_path.parent.resolve()) or identity(path) != artifact_identity:
            raise SystemExit(f"protected artifact mismatch: {label}:{name}")
        artifacts[name] = artifact_identity["sha256"]
    return {"acceptance_sha256": value["sha256"], "run_manifest_sha256": run_identity["sha256"], "artifact_sha256": artifacts}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    run = json.loads((args.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("protected_authority") != EXPECTED_PROTECTED:
        raise SystemExit("protected authority contract mismatch")
    workflow_hash = toolchain(root)
    if run.get("manifest_sha256") != digest(args.manifest) or run.get("toolchain_sha256") != workflow_hash:
        raise SystemExit("input or toolchain identity mismatch")
    verified_source = {}
    for name, value in manifest["files"].items():
        path = (root / value["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file() or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}:
            raise SystemExit(f"source identity mismatch: {name}")
        verified_source[name] = value["sha256"]
    acquisition = json.loads((root / manifest["files"]["acquisition_evidence"]["path"]).read_text(encoding="utf-8"))
    if acquisition.get("raw_read_acquisition", {}).get("status") != "not-performed" or acquisition.get("selected_runs") != SELECTED:
        raise SystemExit("acquisition evidence mismatch")
    reference = json.loads((root / manifest["files"]["reference_manifest"]["path"]).read_text(encoding="utf-8"))["references"]["strain-b"]
    if reference.get("assembly_accession") != "GCA_001746955.1" or reference["annotation_accession"].rsplit(" ", 1)[-1] != manifest["reference_annotation_release"]:
        raise SystemExit("reference manifest mismatch")
    sra = ET.parse(root / manifest["files"]["sra_experiment_packages"]["path"]).getroot()
    ena = {item["run_accession"]: item for item in json.loads((root / manifest["files"]["ena_read_run"]["path"]).read_text(encoding="utf-8"))}
    bioproject = ET.parse(root / manifest["files"]["bioproject"]["path"]).getroot(); archive_id = bioproject.find(".//ArchiveID")
    if archive_id is None or archive_id.get("accession") != "PRJNA1210090" or "engineered k. phaffii" not in " ".join(bioproject.itertext()).lower():
        raise SystemExit("BioProject identity or engineering context mismatch")
    samples = set()
    replicates = set()
    for accession, expected in SELECTED.items():
        package = next(item for item in sra.findall(".//EXPERIMENT_PACKAGE") if item.find("./RUN_SET/RUN") is not None and item.find("./RUN_SET/RUN").get("accession") == accession)
        experiment = package.find("EXPERIMENT"); sample = package.find("SAMPLE"); library = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
        attrs = {item.findtext("TAG", ""): item.findtext("VALUE", "") for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
        biosample = next(item.text for item in sample.findall("./IDENTIFIERS/EXTERNAL_ID") if item.get("namespace") == "BioSample")
        replicate = attrs["replicate"].removeprefix("biological replicate ")
        design = experiment.findtext("./DESIGN/DESIGN_DESCRIPTION", "")
        layout = next(iter(library.find("LIBRARY_LAYOUT"))).tag
        observed = {"experiment": experiment.get("accession"), "sample": biosample, "replicate": replicate}
        if observed != expected or ena[accession]["study_accession"] != "PRJNA1210090" or ena[accession]["secondary_study_accession"] != "SRP557139" or ena[accession]["experiment_accession"] != expected["experiment"] or ena[accession]["sample_accession"] != expected["sample"]:
            raise SystemExit(f"accession mismatch: {accession}")
        if attrs.get("strain") != "GS115" or "wild type" not in design.lower() or "without ferulic acid stress" not in design.lower():
            raise SystemExit(f"identity mismatch: {accession}")
        if (library.findtext("LIBRARY_STRATEGY"), library.findtext("LIBRARY_SOURCE"), library.findtext("LIBRARY_SELECTION"), layout, ena[accession]["instrument_model"]) != ("RNA-Seq", "TRANSCRIPTOMIC", "RANDOM PCR", "PAIRED", "BGISEQ-500"):
            raise SystemExit(f"library mismatch: {accession}")
        if date.fromisoformat(attrs["collection_date"]) <= date.fromisoformat(reference["annotation_accession"].rsplit(" ", 1)[-1]):
            raise SystemExit(f"date mismatch: {accession}")
        samples.add(biosample); replicates.add(replicate)
    if len(samples) != 3 or replicates != {"1", "2", "3"}:
        raise SystemExit("replication mismatch")
    qualification = json.loads((args.run_dir / "metadata_qualification.json").read_text(encoding="utf-8"))
    if qualification.get("gates") != EXPECTED_GATES or qualification.get("scientific_acceptance_status") != "blocked" or qualification.get("eligible_as_independent_transcript_source") is not False:
        raise SystemExit("qualification gate mismatch")
    artifacts = {}
    for name, value in run["artifacts"].items():
        path = args.run_dir / name
        if identity(path) != value:
            raise SystemExit(f"artifact mismatch: {name}")
        artifacts[name] = value["sha256"]
    protected = {name: verify_protected(root, value, name) for name, value in manifest["protected_authority"].items()}
    if protected != run.get("protected_evidence"):
        raise SystemExit("protected evidence mismatch")
    result = {
        "schema_version": 1, "verification_type": "prjna1210090_metadata_qualification", "status": "passed",
        "run_id": run["run_id"], "manifest_sha256": run["manifest_sha256"], "toolchain_sha256": workflow_hash,
        "verified_source_files": verified_source, "verified_artifacts": artifacts,
        "recomputed_gates": EXPECTED_GATES, "scientific_acceptance_status": "blocked",
        "eligible_as_independent_transcript_source": False,
        "protected_evidence": protected,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
