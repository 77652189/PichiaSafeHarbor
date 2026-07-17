from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

from pichia_safe_harbor.independent_transcript_source import toolchain_sha256


REPEATABLE = {
    "source_qualification.json",
    "source_qualification.md",
    "run_manifest.json",
    "verification/independent_check.json",
    "verification/test_evidence.json",
}
EXPECTED = {*REPEATABLE, "verification/execution_receipt.json"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def files(run_dir: Path) -> set[str]:
    return {path.relative_to(run_dir).as_posix() for path in run_dir.rglob("*") if path.is_file()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    if files(args.run_a) != EXPECTED or files(args.run_b) != EXPECTED:
        raise SystemExit("candidate repeatability file set mismatch")
    run_a = json.loads((args.run_a / "run_manifest.json").read_text(encoding="utf-8"))
    run_b = json.loads((args.run_b / "run_manifest.json").read_text(encoding="utf-8"))
    if run_a["run_id"] != run_b["run_id"]:
        raise SystemExit("candidate repeatability run ID mismatch")
    workflow_hash = toolchain_sha256(args.repo_root.resolve())
    if run_a.get("toolchain_sha256") != workflow_hash or run_b.get("toolchain_sha256") != workflow_hash:
        raise SystemExit("candidate repeatability toolchain mismatch")
    hashes = {}
    for relative in sorted(REPEATABLE):
        left = sha256(args.run_a / relative)
        right = sha256(args.run_b / relative)
        if left != right:
            raise SystemExit(f"candidate repeatability mismatch: {relative}")
        hashes[relative] = left
    receipts = {}
    invocation_ids = set()
    coordinator_ids = set()
    run_slots = set()
    for run_dir in (args.run_a, args.run_b):
        receipt_path = run_dir / "verification/execution_receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        steps = receipt.get("workflow_steps", [])
        root = args.repo_root.resolve()
        expected_commands = [
            [sys.executable, str(root / "scripts/run_independent_source_qualification.py"), "--manifest", str(args.manifest.resolve()), "--reference-manifest", str(args.reference_manifest.resolve()), "--output-dir", str(run_dir.resolve()), "--repo-root", str(root)],
            [sys.executable, str(root / "scripts/run_test_evidence.py"), "--output", str((run_dir / "verification/test_evidence.json").resolve())],
            [sys.executable, str(root / "scripts/verify_independent_source.py"), "--run-dir", str(run_dir.resolve()), "--manifest", str(args.manifest.resolve()), "--reference-manifest", str(args.reference_manifest.resolve()), "--repo-root", str(root), "--output", str((run_dir / "verification/independent_check.json").resolve())],
        ]
        commands_valid = len(steps) == 3 and all(step.get("command") == expected for step, expected in zip(steps, expected_commands))
        times_valid = True
        try:
            intervals = [(datetime.fromisoformat(step["started_at_utc"]), datetime.fromisoformat(step["completed_at_utc"])) for step in steps]
            times_valid = all(start <= end for start, end in intervals) and all(intervals[index][1] <= intervals[index + 1][0] for index in range(len(intervals) - 1))
        except (KeyError, TypeError, ValueError):
            times_valid = False
        expected_outputs = {
            relative: {"size_bytes": (run_dir / relative).stat().st_size, "sha256": sha256(run_dir / relative)}
            for relative in sorted(REPEATABLE)
        }
        if (
            receipt.get("evidence_type") != "independent_execution_receipt"
            or receipt.get("candidate_id") != run_a["candidate_id"]
            or receipt.get("run_id") != run_a["run_id"]
            or receipt.get("manifest_sha256") != run_a["manifest_sha256"]
            or receipt.get("reference_manifest_sha256") != run_a["reference_manifest_sha256"]
            or receipt.get("implementation_sha256") != run_a["implementation_sha256"]
            or receipt.get("toolchain_sha256") != workflow_hash
            or not isinstance(receipt.get("command"), list)
            or "run_independent_source_qualification.py" not in " ".join(receipt["command"])
            or receipt.get("workflow_evidence_type") != "orchestrated_independent_source_execution"
            or receipt.get("workflow_status") != "passed"
            or receipt.get("workflow_cwd") != str(args.repo_root.resolve())
            or not commands_valid
            or receipt["command"][1:] != steps[0]["command"][1:]
            or any(step.get("exit_code") != 0 or not step.get("stdout_sha256") for step in steps)
            or not times_valid
            or receipt.get("workflow_outputs") != expected_outputs
        ):
            raise SystemExit(f"candidate execution receipt mismatch: {run_dir.name}")
        invocation_ids.add(receipt.get("invocation_id"))
        coordinator_ids.add(receipt.get("coordinator_id"))
        run_slots.add(receipt.get("run_slot"))
        receipts[run_dir.name] = {"invocation_id": receipt.get("invocation_id"), "identity": {"size_bytes": receipt_path.stat().st_size, "sha256": sha256(receipt_path)}}
    if len(invocation_ids) != 2 or None in invocation_ids or len(coordinator_ids) != 1 or None in coordinator_ids or run_slots != {"run-a", "run-b"}:
        raise SystemExit("candidate repeat runs do not have distinct execution receipts")
    result = {
        "schema_version": 1,
        "verification_type": "independent_transcript_source_repeatability",
        "status": "passed",
        "run_id": run_a["run_id"],
        "toolchain_sha256": workflow_hash,
        "run_directories": [args.run_a.name, args.run_b.name],
        "all_files_identical": True,
        "file_sha256": hashes,
        "execution_receipts": receipts,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
