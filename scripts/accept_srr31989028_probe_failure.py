from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import uuid
from pathlib import Path

import jsonschema


EXPECTED_ARTIFACTS = {
    "alignment.bam", "alignment.bam.bai", "hisat2_summary.txt", "hisat2_stderr.txt", "flagstat.json", "samtools_stats.txt",
    "pair_mapping_summary.json", "splice_junctions.tsv", "coverage_summary.json", "coverage_windows.tsv", "coverage_anomalies.tsv",
    "univec_all_pairs_summary.txt", "univec_both_unmapped_summary.txt", "engineering_impact.json", "capability_report.json", "capability_report.md",
}
EXPECTED_CAPABILITIES = {"coordinate_coverage": "passed", "splice_support": "passed", "strand_support": "unavailable", "boundary_support": "not-authorized", "strain_specific_support": "unavailable", "biological_replication": "single-replicate-only", "unique_multi_unmapped_classification": "passed", "engineering_screen": "failed"}


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): value.update(block)
    return value.hexdigest()


def identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": sha256(path)}


def verified_file(root: Path, value: dict, label: str) -> Path:
    path = (root / value["path"]).resolve()
    if not path.is_relative_to(root) or not path.is_file() or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}: raise SystemExit(f"identity mismatch: {label}")
    return path


def wsl_path(path: Path) -> str:
    value = path.resolve(); return f"/mnt/{value.drive[0].lower()}{value.as_posix()[2:]}"


def semantic_sam_sha256(bam: Path) -> str:
    command = f"set -o pipefail; /usr/bin/samtools quickcheck -v '{wsl_path(bam)}' && /usr/bin/samtools view -h '{wsl_path(bam)}' | sha256sum"
    completed = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], capture_output=True, check=True)
    matches = re.findall(rb"[0-9a-f]{64}", completed.stdout)
    if len(matches) != 1: raise SystemExit("semantic SAM hash missing")
    return matches[0].decode("ascii")


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--repo-root", type=Path, default=Path(".")); parser.add_argument("--manifest", type=Path, required=True); parser.add_argument("--run1", type=Path, required=True); parser.add_argument("--run2", type=Path, required=True); parser.add_argument("--repeatability-output", type=Path, required=True); parser.add_argument("--acceptance-output", type=Path, required=True)
    args = parser.parse_args(); root = args.repo_root.resolve(); manifest_path = args.manifest.resolve(); runs = {"run1": args.run1.resolve(), "run2": args.run2.resolve()}; repeat_path = args.repeatability_output.resolve(); acceptance_path = args.acceptance_output.resolve()
    if any(not path.is_relative_to(root) for path in [manifest_path, *runs.values(), repeat_path, acceptance_path]) or not manifest_path.is_file() or any(not path.is_dir() for path in runs.values()): raise SystemExit("paths must exist inside repo root")
    if repeat_path.exists() or acceptance_path.exists(): raise SystemExit("failure evidence output already exists")
    raw = root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw"
    if any(path.is_file() and path.name.upper().startswith(("SRR31989016", "SRR31989027")) for path in raw.rglob("*")): raise SystemExit("forbidden replicate raw files detected")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")); jsonschema.validate(manifest, json.loads((root / "schemas/srr31989028_probe_manifest.v1.schema.json").read_text(encoding="utf-8")))
    acquisition_path = verified_file(root, manifest["acquisition"]["evidence"], "acquisition evidence"); source_sra = verified_file(root, manifest["acquisition"]["sra"], "source SRA")
    run_values = {name: json.loads((path / "run_manifest.json").read_text(encoding="utf-8")) for name, path in runs.items()}
    if run_values["run1"]["run_id"] != run_values["run2"]["run_id"] or run_values["run1"]["input_manifest_sha256"] != sha256(manifest_path) or run_values["run2"]["input_manifest_sha256"] != sha256(manifest_path): raise SystemExit("run identity mismatch")
    artifacts = {}; capabilities = {}; independents = {}; tests = {}; receipts = {}; invocation_ids = set(); slots = set(); coordinators = set()
    for name, directory in runs.items():
        run = run_values[name]
        if run.get("probe_id") != "prjna1210090-srr31989028-controlled-probe" or run.get("execution_status") != "complete" or run.get("verification_status") != "not_run" or run.get("scientific_acceptance_status") != "blocked" or run.get("parameters") != manifest["parameters"]: raise SystemExit(f"run status/parameter mismatch: {name}")
        if set(run["artifacts"]) != EXPECTED_ARTIFACTS: raise SystemExit(f"artifact set mismatch: {name}")
        artifacts[name] = {}
        for artifact, declared in run["artifacts"].items():
            actual = identity(directory / artifact)
            if actual != declared: raise SystemExit(f"artifact identity mismatch: {name}:{artifact}")
            artifacts[name][artifact] = actual
        capability = json.loads((directory / "capability_report.json").read_text(encoding="utf-8")); engineering = json.loads((directory / "engineering_impact.json").read_text(encoding="utf-8"))
        if capability.get("capabilities") != EXPECTED_CAPABILITIES or engineering.get("widespread_high_depth") is not True or engineering.get("global_univec_signal") is not False: raise SystemExit(f"failure capability mismatch: {name}")
        capabilities[name] = capability["capabilities"]
        independent = json.loads((directory / "verification/independent_check.json").read_text(encoding="utf-8")); test = json.loads((directory / "verification/test_evidence.json").read_text(encoding="utf-8"))
        expected_verified = {artifact: value["sha256"] for artifact, value in artifacts[name].items()}
        expected_stats = {"pair_count": 21957417, "read_count_per_file": 21957417, "total_bases": 6587225100}
        if independent.get("status") != "passed" or independent.get("run_id") != run["run_id"] or independent.get("input_manifest_sha256") != run["input_manifest_sha256"] or independent.get("toolchain_sha256") != run["toolchain_sha256"] or independent.get("verified_artifacts") != expected_verified or independent.get("protected_evidence") != run["protected_evidence"] or independent.get("acquisition_evidence_sha256") != manifest["acquisition"]["evidence"]["sha256"] or independent.get("source_sra_sha256") != manifest["acquisition"]["sra"]["sha256"] or independent.get("recomputed_fastq_statistics") != expected_stats: raise SystemExit(f"independent verification mismatch: {name}")
        if test.get("status") != "passed" or test.get("command") != ["python", "-m", "pytest", "-q"] or test.get("passed_count") != 93 or test.get("toolchain_sha256") != run["toolchain_sha256"]: raise SystemExit(f"test evidence mismatch: {name}")
        independents[name] = independent; tests[name] = test
        receipt_path = directory / "verification/execution_receipt.json"; receipt = json.loads(receipt_path.read_text(encoding="utf-8")); expected_slot = "run-a" if name == "run1" else "run-b"
        try: uuid.UUID(receipt["invocation_id"]); uuid.UUID(receipt["coordinator_id"])
        except (KeyError, TypeError, ValueError): raise SystemExit(f"execution receipt UUID mismatch: {name}")
        if receipt.get("workflow_status") != "passed" or receipt.get("run_slot") != expected_slot or receipt.get("run_id") != run["run_id"] or receipt.get("input_manifest_sha256") != run["input_manifest_sha256"] or receipt.get("toolchain_sha256") != run["toolchain_sha256"] or len(receipt.get("workflow_steps", [])) != 3: raise SystemExit(f"execution receipt mismatch: {name}")
        receipts[name] = receipt; invocation_ids.add(receipt["invocation_id"]); slots.add(receipt["run_slot"]); coordinators.add(receipt["coordinator_id"])
    if run_values["run1"]["toolchain_sha256"] != run_values["run2"]["toolchain_sha256"] or run_values["run1"]["parameters"] != run_values["run2"]["parameters"] or run_values["run1"]["protected_evidence"] != run_values["run2"]["protected_evidence"] or len(invocation_ids) != 2 or slots != {"run-a", "run-b"} or len(coordinators) != 1: raise SystemExit("independent execution binding mismatch")
    identical = {}; different = {}
    for artifact in sorted(EXPECTED_ARTIFACTS):
        left = artifacts["run1"][artifact]["sha256"]; right = artifacts["run2"][artifact]["sha256"]
        (identical if left == right else different)[artifact] = left if left == right else {"run1": left, "run2": right}
    if set(different) != {"alignment.bam", "alignment.bam.bai"}: raise SystemExit(f"unexpected repeatability difference: {sorted(different)}")
    semantic = {name: semantic_sam_sha256(directory / "alignment.bam") for name, directory in runs.items()}
    semantic_identical = len(set(semantic.values())) == 1
    repeatability = {"schema_version": 1, "verification_type": "srr31989028_controlled_probe_repeatability", "status": "failed", "run_id": run_values["run1"]["run_id"], "compared_artifact_count": len(EXPECTED_ARTIFACTS), "identical_artifact_sha256": identical, "different_artifact_sha256": different, "sam_record_stream_sha256": semantic, "sam_record_stream_identical": semantic_identical, "byte_identical_main_artifacts": False, "execution_receipts": {name: {"identity": identity(runs[name] / "verification/execution_receipt.json"), "invocation_id": receipts[name]["invocation_id"], "run_slot": receipts[name]["run_slot"], "coordinator_id": receipts[name]["coordinator_id"]} for name in runs}, "failure_reason": "BAM record stream ordering/bytes and derived BAI offsets differ across independent runs"}
    repeat_schema = json.loads((root / "schemas/srr31989028_probe_failure_repeatability.v1.schema.json").read_text(encoding="utf-8")); jsonschema.validate(repeatability, repeat_schema)
    repeat_path.parent.mkdir(parents=True, exist_ok=True); repeat_temp = repeat_path.with_name(repeat_path.name + ".tmp"); acceptance_temp = acceptance_path.with_name(acceptance_path.name + ".tmp"); repeat_temp.write_text(json.dumps(repeatability, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run1_tests = identity(runs["run1"] / "verification/test_evidence.json") | {"passed_count": tests["run1"]["passed_count"]}; run2_tests = identity(runs["run2"] / "verification/test_evidence.json") | {"passed_count": tests["run2"]["passed_count"]}
    acceptance = {"schema_version": 1, "acceptance_version": "srr31989028-controlled-probe-failure-v1", "run_id": run_values["run1"]["run_id"], "probe_id": "prjna1210090-srr31989028-controlled-probe", "probe_acceptance_status": "failed", "execution_status": "complete", "verification_status": "passed", "scientific_acceptance_status": "blocked", "failure_reasons": ["widespread-high-depth-engineering-screen", "bam-bai-byte-repeatability"], "capabilities": EXPECTED_CAPABILITIES, "strand_support": "unavailable", "single_replicate_only": True, "three_replicate_evidence": False, "formal_boundary_track": False, "input_manifest": identity(manifest_path), "acquisition_evidence": identity(acquisition_path), "source_sra": identity(source_sra), "run_manifests": {name: identity(directory / "run_manifest.json") for name, directory in runs.items()}, "artifacts_by_run": artifacts, "verification_evidence": {"run1_independent": identity(runs["run1"] / "verification/independent_check.json"), "run2_independent": identity(runs["run2"] / "verification/independent_check.json"), "run1_tests": run1_tests, "run2_tests": run2_tests}, "execution_receipts": {name: identity(runs[name] / "verification/execution_receipt.json") for name in runs}, "repeatability_evidence": identity(repeat_temp), "protected_evidence": run_values["run1"]["protected_evidence"], "acceptance_tool_sha256": sha256(Path(__file__)), "stop_line": "single SRR31989028 probe failed; remain blocked and do not download other replicates, create a boundary track, generate candidates or thresholds, or enter Slice 1"}
    schema = json.loads((root / "schemas/srr31989028_probe_failure_acceptance.v1.schema.json").read_text(encoding="utf-8")); jsonschema.validate(acceptance, schema)
    acceptance_path.parent.mkdir(parents=True, exist_ok=True); acceptance_temp.write_text(json.dumps(acceptance, indent=2, sort_keys=True) + "\n", encoding="utf-8"); os.replace(repeat_temp, repeat_path); os.replace(acceptance_temp, acceptance_path)
    print(json.dumps({"run_id": acceptance["run_id"], "probe_acceptance_status": "failed", "statuses": ["complete", "passed", "blocked"], "failure_reasons": acceptance["failure_reasons"], "sam_record_stream_sha256": semantic}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
