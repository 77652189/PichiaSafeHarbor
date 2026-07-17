from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

import jsonschema

from run_srr31989028_controlled_probe import EXPECTED_ARTIFACTS, identity, sha256, toolchain_sha256


REPEATABLE = {*EXPECTED_ARTIFACTS, "run_manifest.json", "verification/independent_check.json", "verification/test_evidence.json"}
PREACCEPTANCE = {*REPEATABLE, "verification/execution_receipt.json"}
REQUIRED_CAPABILITIES = {
    "coordinate_coverage": "passed",
    "splice_support": "passed",
    "strand_support": "unavailable",
    "boundary_support": "not-authorized",
    "strain_specific_support": "unavailable",
    "biological_replication": "single-replicate-only",
    "unique_multi_unmapped_classification": "passed",
    "engineering_screen": "passed-with-explicit-local-exclusions",
}


def run_step(command: list[str], root: Path) -> tuple[dict, bytes, bytes]:
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(command, cwd=root, capture_output=True, check=False)
    ended = datetime.now(timezone.utc).isoformat()
    evidence = {"command": command, "cwd": str(root), "started_at_utc": started, "completed_at_utc": ended, "exit_code": completed.returncode, "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(), "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest()}
    if completed.returncode != 0:
        sys.stdout.buffer.write(completed.stdout); sys.stderr.buffer.write(completed.stderr)
        raise SystemExit(f"workflow step failed: {command}")
    return evidence, completed.stdout, completed.stderr


def file_set(path: Path) -> set[str]:
    return {item.relative_to(path).as_posix() for item in path.rglob("*") if item.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True); parser.add_argument("--run-a", type=Path, required=True); parser.add_argument("--run-b", type=Path, required=True); parser.add_argument("--repeatability-output", type=Path, required=True); parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args(); root = args.repo_root.resolve(); manifest_path = args.manifest.resolve(); run_a_path = args.run_a.resolve(); run_b_path = args.run_b.resolve(); repeat_path = args.repeatability_output.resolve()
    if not manifest_path.is_relative_to(root) or not manifest_path.is_file(): raise SystemExit("manifest must be inside repo root")
    if any(not path.is_relative_to(root) for path in (run_a_path, run_b_path, repeat_path)): raise SystemExit("workflow outputs must be inside repo root")
    if len({run_a_path, run_b_path, repeat_path}) != 3 or any(path.exists() for path in (run_a_path, run_b_path, repeat_path)): raise SystemExit("workflow outputs must be distinct and absent")
    manifest_schema = json.loads((root / "schemas/srr31989028_probe_manifest.v1.schema.json").read_text(encoding="utf-8")); jsonschema.validate(json.loads(manifest_path.read_text(encoding="utf-8")), manifest_schema)
    workflow_hash = toolchain_sha256(root); coordinator_id = str(uuid.uuid4())
    for slot, run_dir in (("run-a", run_a_path), ("run-b", run_b_path)):
        producer = [sys.executable, str(root / "scripts/run_srr31989028_controlled_probe.py"), "--manifest", str(manifest_path), "--output-dir", str(run_dir), "--repo-root", str(root)]
        tests = [sys.executable, "-m", "pytest", "-q"]
        verifier = [sys.executable, str(root / "scripts/verify_srr31989028_controlled_probe.py"), "--run-dir", str(run_dir), "--manifest", str(manifest_path), "--repo-root", str(root), "--output", str(run_dir / "verification/independent_check.json")]
        step1, _, _ = run_step(producer, root); step2, stdout, stderr = run_step(tests, root)
        match = re.search(rb"(\d+) passed", stdout + stderr)
        if not match: raise SystemExit("pytest pass count missing")
        test_evidence = {"schema_version": 1, "evidence_type": "automated_tests", "status": "passed", "command": ["python", "-m", "pytest", "-q"], "passed_count": int(match.group(1)), "toolchain_sha256": workflow_hash}
        test_path = run_dir / "verification/test_evidence.json"; test_path.write_text(json.dumps(test_evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        step3, _, _ = run_step(verifier, root)
        receipt_path = run_dir / "verification/execution_receipt.json"; receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt.update({"workflow_evidence_type": "orchestrated_srr31989028_controlled_probe", "workflow_status": "passed", "workflow_cwd": str(root), "coordinator_id": coordinator_id, "run_slot": slot, "workflow_completed_at_utc": datetime.now(timezone.utc).isoformat(), "workflow_steps": [step1, step2, step3], "workflow_outputs": {name: identity(run_dir / name) for name in sorted(REPEATABLE)}})
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_a = json.loads((run_a_path / "run_manifest.json").read_text(encoding="utf-8")); run_b = json.loads((run_b_path / "run_manifest.json").read_text(encoding="utf-8"))
    if run_a != run_b or run_a.get("toolchain_sha256") != workflow_hash: raise SystemExit("repeat run manifest mismatch")
    hashes = {}
    for name in sorted(REPEATABLE):
        left = sha256(run_a_path / name); right = sha256(run_b_path / name)
        if left != right: raise SystemExit(f"repeatability mismatch: {name}")
        hashes[name] = left
    receipts = {}; invocation_ids = set(); slots = set(); coordinators = set()
    for run_dir in (run_a_path, run_b_path):
        receipt_path = run_dir / "verification/execution_receipt.json"; receipt = json.loads(receipt_path.read_text(encoding="utf-8")); steps = receipt.get("workflow_steps", [])
        expected = [[sys.executable, str(root / "scripts/run_srr31989028_controlled_probe.py"), "--manifest", str(manifest_path), "--output-dir", str(run_dir), "--repo-root", str(root)], [sys.executable, "-m", "pytest", "-q"], [sys.executable, str(root / "scripts/verify_srr31989028_controlled_probe.py"), "--run-dir", str(run_dir), "--manifest", str(manifest_path), "--repo-root", str(root), "--output", str(run_dir / "verification/independent_check.json")]]
        if receipt.get("workflow_evidence_type") != "orchestrated_srr31989028_controlled_probe" or receipt.get("workflow_status") != "passed" or receipt.get("workflow_cwd") != str(root) or receipt.get("run_id") != run_a["run_id"] or receipt.get("input_manifest_sha256") != run_a["input_manifest_sha256"] or receipt.get("toolchain_sha256") != workflow_hash or len(steps) != 3 or any(step.get("command") != command or step.get("cwd") != str(root) or step.get("exit_code") != 0 for step, command in zip(steps, expected)) or receipt.get("workflow_outputs") != {name: identity(run_dir / name) for name in sorted(REPEATABLE)}: raise SystemExit(f"execution receipt mismatch: {run_dir.name}")
        try: uuid.UUID(receipt["invocation_id"]); uuid.UUID(receipt["coordinator_id"])
        except (KeyError, TypeError, ValueError): raise SystemExit("execution receipt UUID mismatch")
        invocation_ids.add(receipt["invocation_id"]); slots.add(receipt["run_slot"]); coordinators.add(receipt["coordinator_id"]); receipts[run_dir.name] = {"invocation_id": receipt["invocation_id"], "identity": identity(receipt_path)}
    if len(invocation_ids) != 2 or slots != {"run-a", "run-b"} or coordinators != {coordinator_id}: raise SystemExit("repeat execution independence mismatch")
    repeatability = {"schema_version": 1, "verification_type": "srr31989028_controlled_probe_repeatability", "status": "passed", "run_id": run_a["run_id"], "toolchain_sha256": workflow_hash, "run_directories": [run_a_path.name, run_b_path.name], "compared_file_count": len(hashes), "all_files_identical": True, "file_sha256": hashes, "execution_receipts": receipts}
    repeat_path.parent.mkdir(parents=True, exist_ok=True); repeat_path.write_text(json.dumps(repeatability, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if file_set(run_a_path) != PREACCEPTANCE or file_set(run_b_path) != PREACCEPTANCE: raise SystemExit("preacceptance file set mismatch")
    capability = json.loads((run_a_path / "capability_report.json").read_text(encoding="utf-8")); independent = json.loads((run_a_path / "verification/independent_check.json").read_text(encoding="utf-8")); tests = json.loads((run_a_path / "verification/test_evidence.json").read_text(encoding="utf-8"))
    if capability.get("capabilities") != REQUIRED_CAPABILITIES or capability.get("strand_direction") != "unavailable" or capability.get("single_replicate_only") is not True or capability.get("three_replicate_evidence") is not False or capability.get("formal_boundary_track") is not False: raise SystemExit("capability acceptance gates failed")
    if independent.get("status") != "passed" or independent.get("verified_artifacts") != {name: value["sha256"] for name, value in run_a["artifacts"].items()} or independent.get("protected_evidence") != run_a.get("protected_evidence"): raise SystemExit("independent verification gate failed")
    if tests.get("status") != "passed" or tests.get("passed_count", 0) < 1 or tests.get("toolchain_sha256") != workflow_hash: raise SystemExit("automated test gate failed")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    raw = root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw"
    if any(path.is_file() and path.name.upper().startswith(("SRR31989016", "SRR31989027")) for path in raw.rglob("*")): raise SystemExit("forbidden replicate raw files detected")
    current_acquisition = identity(root / manifest["acquisition"]["evidence"]["path"]); current_sra = identity(root / manifest["acquisition"]["sra"]["path"])
    if current_acquisition != {key: manifest["acquisition"]["evidence"][key] for key in ("size_bytes", "sha256")} or current_sra != {key: manifest["acquisition"]["sra"][key] for key in ("size_bytes", "sha256")} or independent.get("acquisition_evidence_sha256") != current_acquisition["sha256"] or independent.get("source_sra_sha256") != current_sra["sha256"] or independent.get("recomputed_fastq_statistics") != {"pair_count": 21957417, "read_count_per_file": 21957417, "total_bases": 6587225100}: raise SystemExit("independent acquisition evidence mismatch")
    acceptance = {"schema_version": 1, "acceptance_version": "srr31989028-controlled-probe-v1", "run_id": run_a["run_id"], "probe_id": run_a["probe_id"], "probe_acceptance_status": "passed", "execution_status": "complete", "verification_status": "passed", "scientific_acceptance_status": "blocked", "capabilities": REQUIRED_CAPABILITIES, "strand_support": "unavailable", "single_replicate_only": True, "three_replicate_evidence": False, "formal_boundary_track": False, "input_manifest": identity(manifest_path), "acquisition_evidence": identity(root / manifest["acquisition"]["evidence"]["path"]), "source_sra": identity(root / manifest["acquisition"]["sra"]["path"]), "run_manifest": identity(run_a_path / "run_manifest.json"), "artifacts": run_a["artifacts"], "verification_evidence": {"independent": identity(run_a_path / "verification/independent_check.json"), "automated_tests": identity(run_a_path / "verification/test_evidence.json") | {"passed_count": tests["passed_count"]}, "repeatability": identity(repeat_path)}, "protected_evidence": run_a["protected_evidence"], "stop_line": "single SRR31989028 probe complete; remain blocked and do not download other replicates, create a boundary track, generate candidates or thresholds, or enter Slice 1"}
    schema = json.loads((root / "schemas/srr31989028_probe_acceptance.v1.schema.json").read_text(encoding="utf-8")); jsonschema.validate(acceptance, schema)
    (run_a_path / "acceptance_manifest.json").write_text(json.dumps(acceptance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"run_id": run_a["run_id"], "probe_acceptance_status": "passed", "statuses": ["complete", "passed", "blocked"], "test_count": tests["passed_count"], "repeatable_files": len(hashes)}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
