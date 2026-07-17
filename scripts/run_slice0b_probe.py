from __future__ import annotations

import argparse
import json
from pathlib import Path

from pichia_safe_harbor.transcript_probe import run_transcript_probe


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--toolchain", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    result = run_transcript_probe(args.sources, args.toolchain, args.reference_manifest, args.output_dir, args.repo_root)
    print(json.dumps({"run_id": result["run_id"], "execution_status": result["execution_status"], "scientific_acceptance_status": result["scientific_acceptance_status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
