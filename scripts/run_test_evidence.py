from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

from pichia_safe_harbor.pipeline import _implementation_hash
from pichia_safe_harbor.independent_transcript_source import toolchain_sha256


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    command = [sys.executable, "-m", "pytest", "-q"]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    output = completed.stdout + completed.stderr
    match = re.search(r"(\d+) passed", output)
    result = {
        "schema_version": 1,
        "evidence_type": "automated_tests",
        "status": "passed" if completed.returncode == 0 and match else "failed",
        "command": ["python", "-m", "pytest", "-q"],
        "passed_count": int(match.group(1)) if match else 0,
        "implementation_sha256": _implementation_hash(),
        "toolchain_sha256": toolchain_sha256(Path(__file__).resolve().parents[1]),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output, end="")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
