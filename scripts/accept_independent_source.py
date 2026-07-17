from __future__ import annotations

import argparse
import json
from pathlib import Path

from pichia_safe_harbor.independent_transcript_source import create_independent_source_acceptance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--peer-run-dir", type=Path, required=True)
    parser.add_argument("--independent-check", type=Path, required=True)
    parser.add_argument("--test-evidence", type=Path, required=True)
    parser.add_argument("--repeatability-evidence", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    args = parser.parse_args()
    result = create_independent_source_acceptance(args.run_dir, args.peer_run_dir, args.independent_check, args.test_evidence, args.repeatability_evidence, args.repo_root, args.manifest, args.reference_manifest)
    print(json.dumps({key: result[key] for key in ("run_id", "execution_status", "verification_status", "scientific_acceptance_status")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
