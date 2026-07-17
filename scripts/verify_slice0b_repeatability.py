from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_FILES = {
    "run_manifest.json",
    "slice0b_recommendation.md",
    "source_applicability_matrix.json",
    "source_applicability_matrix.tsv",
    "transcript_evidence_probe.json",
    "transcript_evidence_probe.tsv",
    "transcript_quality_report.md",
    "verification/independent_check.json",
    "verification/test_evidence.json",
}


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
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    if files(args.run_a) != EXPECTED_FILES or files(args.run_b) != EXPECTED_FILES:
        raise SystemExit("Slice 0B repeatability file set mismatch")
    manifest_a = json.loads((args.run_a / "run_manifest.json").read_text(encoding="utf-8"))
    manifest_b = json.loads((args.run_b / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest_a["run_id"] != manifest_b["run_id"]:
        raise SystemExit("Slice 0B repeatability run ID mismatch")
    hashes = {}
    for relative in sorted(EXPECTED_FILES):
        left = sha256(args.run_a / Path(relative))
        right = sha256(args.run_b / Path(relative))
        if left != right:
            raise SystemExit(f"Slice 0B repeatability mismatch: {relative}")
        hashes[relative] = left
    result = {
        "schema_version": 1,
        "verification_type": "slice0b_repeatability",
        "status": "passed",
        "run_id": manifest_a["run_id"],
        "run_directories": [args.run_a.name, args.run_b.name],
        "compared_file_count": len(hashes),
        "all_files_identical": True,
        "file_sha256": hashes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
