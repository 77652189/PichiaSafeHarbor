from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

from pichia_safe_harbor.independent_transcript_source import toolchain_sha256


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _run(command: list[str], cwd: Path) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(command, cwd=cwd, capture_output=True, check=False)
    ended = datetime.now(timezone.utc).isoformat()
    if completed.returncode != 0:
        sys.stdout.buffer.write(completed.stdout)
        sys.stderr.buffer.write(completed.stderr)
        raise SystemExit(f"workflow step failed: {command[1]}")
    return {
        "command": command,
        "started_at_utc": started,
        "completed_at_utc": ended,
        "exit_code": completed.returncode,
        "stdout_sha256": _sha256_bytes(completed.stdout),
        "stderr_sha256": _sha256_bytes(completed.stderr),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument("--repeatability-output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    outputs = [args.run_a, args.run_b, args.repeatability_output]
    if any(path.exists() for path in outputs):
        raise SystemExit("workflow outputs must not already exist")
    coordinator_id = str(uuid.uuid4())
    workflow_hash = toolchain_sha256(root)
    for slot, run_dir in (("run-a", args.run_a), ("run-b", args.run_b)):
        commands = [
            [sys.executable, str(root / "scripts/run_independent_source_qualification.py"), "--manifest", str(args.manifest.resolve()), "--reference-manifest", str(args.reference_manifest.resolve()), "--output-dir", str(run_dir.resolve()), "--repo-root", str(root)],
            [sys.executable, str(root / "scripts/run_test_evidence.py"), "--output", str((run_dir / "verification/test_evidence.json").resolve())],
            [sys.executable, str(root / "scripts/verify_independent_source.py"), "--run-dir", str(run_dir.resolve()), "--manifest", str(args.manifest.resolve()), "--reference-manifest", str(args.reference_manifest.resolve()), "--repo-root", str(root), "--output", str((run_dir / "verification/independent_check.json").resolve())],
        ]
        steps = [_run(command, root) for command in commands]
        receipt_path = run_dir / "verification/execution_receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt.update({
            "workflow_evidence_type": "orchestrated_independent_source_execution",
            "workflow_status": "passed",
            "workflow_cwd": str(root),
            "coordinator_id": coordinator_id,
            "run_slot": slot,
            "toolchain_sha256": workflow_hash,
            "workflow_steps": steps,
            "workflow_outputs": {
                relative: _identity(run_dir / relative)
                for relative in (
                    "run_manifest.json", "source_qualification.json", "source_qualification.md",
                    "verification/independent_check.json", "verification/test_evidence.json",
                )
            },
        })
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    repeat_command = [sys.executable, str(root / "scripts/verify_independent_source_repeatability.py"), "--run-a", str(args.run_a.resolve()), "--run-b", str(args.run_b.resolve()), "--output", str(args.repeatability_output.resolve()), "--repo-root", str(root), "--manifest", str(args.manifest.resolve()), "--reference-manifest", str(args.reference_manifest.resolve())]
    _run(repeat_command, root)
    print(json.dumps({"coordinator_id": coordinator_id, "run_a": str(args.run_a), "run_b": str(args.run_b), "repeatability": str(args.repeatability_output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
