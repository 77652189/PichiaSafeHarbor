from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


SOURCE_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR31989028/SRR31989028"
SRA_SIZE = 3_931_837_205
SRA_MD5 = "8c8a58880254746890e10c68b402c875"


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def file_identity(path: Path, root: Path) -> dict:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise ValueError(f"file must be inside repo root: {path}")
    return {"path": resolved.relative_to(root).as_posix(), "size_bytes": resolved.stat().st_size, "sha256": digest(resolved, "sha256")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--sra", type=Path, required=True)
    parser.add_argument("--aria2c", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    raw = (root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw").resolve()
    sra = args.sra.resolve()
    output = args.output.resolve()
    aria2c = args.aria2c.resolve()
    if not sra.is_relative_to(raw) or not output.is_relative_to(raw) or output.exists() or sra.with_suffix(sra.suffix + ".aria2").exists():
        raise SystemExit("download is incomplete or receipt path is invalid")
    actual = {"size_bytes": sra.stat().st_size, "md5": digest(sra, "md5"), "sha256": digest(sra, "sha256")}
    if actual["size_bytes"] != SRA_SIZE or actual["md5"] != SRA_MD5:
        raise SystemExit("downloaded normalized SRA identity mismatch")
    version = subprocess.run([str(aria2c), "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    argv = [
        str(aria2c), "--continue=true", "--max-connection-per-server=8", "--split=8", "--min-split-size=16M",
        "--file-allocation=none", "--auto-file-renaming=false", "--allow-overwrite=false", "--summary-interval=15",
        "--console-log-level=notice", f"--dir={raw}", "--out=SRR31989028.sra", SOURCE_URL,
    ]
    receipt = {
        "schema_version": 1,
        "download_status": "complete",
        "source_url": SOURCE_URL,
        "declared_size_bytes": SRA_SIZE,
        "declared_md5": SRA_MD5,
        "actual_sra": actual,
        "tool": file_identity(aria2c, root),
        "tool_version": version,
        "argv": argv,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"download_status": "complete", "actual_sra": actual}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
