from __future__ import annotations

import argparse
import json
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

from pichia_safe_harbor.independent_transcript_source import qualify_independent_transcript_source
from pichia_safe_harbor.independent_transcript_source import toolchain_sha256


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    result = qualify_independent_transcript_source(args.manifest, args.reference_manifest, args.output_dir, args.repo_root)
    receipt_path = args.output_dir / "verification/execution_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps({
        "schema_version": 1,
        "evidence_type": "independent_execution_receipt",
        "candidate_id": "prjna604658-srr11011828",
        "run_id": result["run_id"],
        "invocation_id": str(uuid.uuid4()),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": ["python", *sys.argv],
        "manifest_sha256": result["manifest_sha256"],
        "reference_manifest_sha256": result["reference_manifest_sha256"],
        "implementation_sha256": result["implementation_sha256"],
        "toolchain_sha256": toolchain_sha256(Path(__file__).resolve().parents[1]),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"run_id": result["run_id"], "execution_status": result["execution_status"], "scientific_acceptance_status": result["scientific_acceptance_status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
