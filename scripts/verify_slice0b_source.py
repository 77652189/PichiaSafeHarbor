from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from pichia_safe_harbor.pipeline import _implementation_hash
from pichia_safe_harbor.transcript_qualification import (
    load_transcript_source_manifest,
    qualify_transcript_source,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    manifest = load_transcript_source_manifest(args.manifest)
    sources = [qualify_transcript_source(source, args.repo_root) for source in manifest["sources"]]
    result = {
        "schema_version": 1,
        "verification_type": "slice0b_transcript_source_identity",
        "status": "passed",
        "primary_coordinate_space": manifest["primary_coordinate_space"],
        "exact_target_strain_coordinates": manifest["exact_target_strain_coordinates"],
        "manifest_sha256": sha256(args.manifest),
        "acquisition_evidence_sha256": sha256(args.repo_root / manifest["acquisition_evidence"]),
        "implementation_sha256": _implementation_hash(),
        "sources": sources,
        "slice0b_gate_status": "mapping-pending",
        "scientific_acceptance_status": "blocked",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
