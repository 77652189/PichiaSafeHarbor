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

from run_prjna1210090_metadata import CANDIDATE_ID, GATES, SELECTED, identity, sha256, toolchain_sha256


REPEATABLE = {"metadata_qualification.json", "metadata_qualification.md", "run_manifest.json", "verification/independent_check.json", "verification/test_evidence.json"}
PREACCEPTANCE = {*REPEATABLE, "verification/execution_receipt.json"}


def run_step(command: list[str], root: Path) -> tuple[dict, bytes, bytes]:
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(command, cwd=root, capture_output=True, check=False)
    ended = datetime.now(timezone.utc).isoformat()
    evidence = {
        "command": command, "started_at_utc": started, "completed_at_utc": ended,
        "exit_code": completed.returncode,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
    }
    if completed.returncode != 0:
        sys.stdout.buffer.write(completed.stdout); sys.stderr.buffer.write(completed.stderr)
        raise SystemExit(f"workflow step failed: {command[1]}")
    return evidence, completed.stdout, completed.stderr


def file_set(path: Path) -> set[str]:
    return {item.relative_to(path).as_posix() for item in path.rglob("*") if item.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument("--repeatability-output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve(); manifest_path = args.manifest.resolve()
    if any(path.exists() for path in (args.run_a, args.run_b, args.repeatability_output)):
        raise SystemExit("workflow outputs already exist")
    workflow_hash = toolchain_sha256(root)
    coordinator_id = str(uuid.uuid4())
    for slot, run_dir in (("run-a", args.run_a), ("run-b", args.run_b)):
        producer = [sys.executable, str(root / "scripts/run_prjna1210090_metadata.py"), "--manifest", str(manifest_path), "--output-dir", str(run_dir.resolve()), "--repo-root", str(root)]
        verifier = [sys.executable, str(root / "scripts/verify_prjna1210090_metadata.py"), "--run-dir", str(run_dir.resolve()), "--manifest", str(manifest_path), "--repo-root", str(root), "--output", str((run_dir / "verification/independent_check.json").resolve())]
        step1, _, _ = run_step(producer, root)
        test_command = [sys.executable, "-m", "pytest", "-q"]
        step2, test_stdout, test_stderr = run_step(test_command, root)
        match = re.search(rb"(\d+) passed", test_stdout + test_stderr)
        if not match:
            raise SystemExit("pytest pass count missing")
        test_evidence = {
            "schema_version": 1, "evidence_type": "automated_tests", "status": "passed",
            "command": ["python", "-m", "pytest", "-q"], "passed_count": int(match.group(1)),
            "toolchain_sha256": workflow_hash,
        }
        test_path = run_dir / "verification/test_evidence.json"
        test_path.write_text(json.dumps(test_evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        step3, _, _ = run_step(verifier, root)
        receipt_path = run_dir / "verification/execution_receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt.update({
            "workflow_evidence_type": "orchestrated_prjna1210090_metadata_execution", "workflow_status": "passed",
            "workflow_cwd": str(root), "coordinator_id": coordinator_id, "run_slot": slot,
            "workflow_completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "workflow_steps": [step1, step2, step3],
            "workflow_outputs": {relative: identity(run_dir / relative) for relative in sorted(REPEATABLE)},
        })
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_a = json.loads((args.run_a / "run_manifest.json").read_text(encoding="utf-8")); run_b = json.loads((args.run_b / "run_manifest.json").read_text(encoding="utf-8"))
    if run_a != run_b or run_a["toolchain_sha256"] != workflow_hash:
        raise SystemExit("repeat run manifest mismatch")
    hashes = {}
    for relative in sorted(REPEATABLE):
        left = sha256(args.run_a / relative); right = sha256(args.run_b / relative)
        if left != right:
            raise SystemExit(f"repeatability mismatch: {relative}")
        hashes[relative] = left
    receipts = {}
    invocation_ids = set(); slots = set(); coordinators = set()
    for run_dir in (args.run_a, args.run_b):
        path = run_dir / "verification/execution_receipt.json"; value = json.loads(path.read_text(encoding="utf-8"))
        expected_commands = [
            [sys.executable, str(root / "scripts/run_prjna1210090_metadata.py"), "--manifest", str(manifest_path), "--output-dir", str(run_dir.resolve()), "--repo-root", str(root)],
            [sys.executable, "-m", "pytest", "-q"],
            [sys.executable, str(root / "scripts/verify_prjna1210090_metadata.py"), "--run-dir", str(run_dir.resolve()), "--manifest", str(manifest_path), "--repo-root", str(root), "--output", str((run_dir / "verification/independent_check.json").resolve())],
        ]
        steps = value.get("workflow_steps", [])
        try:
            intervals = [(datetime.fromisoformat(step["started_at_utc"]), datetime.fromisoformat(step["completed_at_utc"])) for step in steps]
            times_valid = len(intervals) == 3 and all(start <= end for start, end in intervals) and all(intervals[index][1] <= intervals[index + 1][0] for index in range(2))
            producer_completed = datetime.fromisoformat(value["producer_completed_at_utc"])
            workflow_completed = datetime.fromisoformat(value["workflow_completed_at_utc"])
            times_valid = times_valid and intervals[0][0] <= producer_completed <= intervals[0][1] and workflow_completed >= intervals[-1][1]
        except (KeyError, TypeError, ValueError):
            times_valid = False
        if (
            value.get("schema_version") != 1 or value.get("evidence_type") != "independent_execution_receipt"
            or value.get("candidate_id") != CANDIDATE_ID or value.get("run_id") != run_a["run_id"]
            or value.get("workflow_evidence_type") != "orchestrated_prjna1210090_metadata_execution"
            or value.get("workflow_status") != "passed" or value.get("workflow_cwd") != str(root)
            or value.get("manifest_sha256") != run_a["manifest_sha256"] or value.get("toolchain_sha256") != workflow_hash
            or len(steps) != 3 or any(step.get("command") != expected for step, expected in zip(steps, expected_commands))
            or any(step.get("exit_code") != 0 or not step.get("stdout_sha256") for step in steps) or not times_valid
            or value.get("workflow_outputs") != {relative: identity(run_dir / relative) for relative in sorted(REPEATABLE)}
        ):
            raise SystemExit(f"execution receipt content mismatch: {run_dir.name}")
        try:
            uuid.UUID(value["invocation_id"]); uuid.UUID(value["coordinator_id"])
        except (KeyError, ValueError, TypeError):
            raise SystemExit(f"execution receipt UUID mismatch: {run_dir.name}")
        invocation_ids.add(value["invocation_id"]); slots.add(value["run_slot"]); coordinators.add(value["coordinator_id"])
        receipts[run_dir.name] = {"invocation_id": value["invocation_id"], "identity": identity(path)}
    if len(invocation_ids) != 2 or slots != {"run-a", "run-b"} or coordinators != {coordinator_id}:
        raise SystemExit("execution receipt independence mismatch")
    repeatability = {
        "schema_version": 1, "verification_type": "prjna1210090_metadata_repeatability", "status": "passed",
        "run_id": run_a["run_id"], "toolchain_sha256": workflow_hash, "run_directories": [args.run_a.name, args.run_b.name],
        "all_files_identical": True, "file_sha256": hashes, "execution_receipts": receipts,
    }
    args.repeatability_output.parent.mkdir(parents=True, exist_ok=True)
    args.repeatability_output.write_text(json.dumps(repeatability, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if file_set(args.run_a) != PREACCEPTANCE or file_set(args.run_b) != PREACCEPTANCE:
        raise SystemExit("preacceptance file set mismatch")
    qualification = json.loads((args.run_a / "metadata_qualification.json").read_text(encoding="utf-8"))
    independent = json.loads((args.run_a / "verification/independent_check.json").read_text(encoding="utf-8"))
    tests = json.loads((args.run_a / "verification/test_evidence.json").read_text(encoding="utf-8"))
    if qualification["gates"] != GATES or qualification["scientific_acceptance_status"] != "blocked" or qualification["eligible_as_independent_transcript_source"] is not False:
        raise SystemExit("qualification overclaim")
    if independent.get("status") != "passed" or independent.get("recomputed_gates") != GATES or independent.get("toolchain_sha256") != workflow_hash:
        raise SystemExit("independent verification mismatch")
    if independent.get("verified_source_files") != {name: value["sha256"] for name, value in run_a["source_files"].items()} or independent.get("protected_evidence") != run_a.get("protected_evidence"):
        raise SystemExit("independent source or protected evidence mismatch")
    if tests.get("status") != "passed" or tests.get("passed_count", 0) < 70 or tests.get("toolchain_sha256") != workflow_hash:
        raise SystemExit("test evidence mismatch")
    acceptance = {
        "schema_version": 1, "acceptance_version": "prjna1210090-metadata-v1", "run_id": run_a["run_id"],
        "candidate_id": CANDIDATE_ID, "execution_status": "complete", "verification_status": "passed",
        "scientific_acceptance_status": "blocked", "eligible_as_independent_transcript_source": False,
        "toolchain_sha256": workflow_hash, "qualification_gates": GATES,
        "blocking_gates": {key: value for key, value in GATES.items() if value != "passed"},
        "run_manifest": identity(args.run_a / "run_manifest.json"), "artifacts": run_a["artifacts"],
        "verification_evidence": {
            "independent": identity(args.run_a / "verification/independent_check.json"),
            "automated_tests": identity(args.run_a / "verification/test_evidence.json") | {"passed_count": tests["passed_count"]},
            "repeatability": identity(args.repeatability_output),
        },
        "stop_line": "remain blocked; do not download raw reads, run coordinate mapping, reannotate, or enter Slice 1",
    }
    acceptance_schema = json.loads((root / "schemas/prjna1210090_metadata_acceptance.v1.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(acceptance, acceptance_schema)
    (args.run_a / "acceptance_manifest.json").write_text(json.dumps(acceptance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"run_id": run_a["run_id"], "statuses": ["complete", "passed", "blocked"], "test_count": tests["passed_count"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
