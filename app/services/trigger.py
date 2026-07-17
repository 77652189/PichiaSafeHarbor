from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pichia_safe_harbor.errors import ContractError  # noqa: E402
from pichia_safe_harbor.pipeline import run_baseline  # noqa: E402
from pichia_safe_harbor.reference import (  # noqa: E402
    load_reference_manifest,
    reference_by_key,
    verify_reference_files,
)
from pichia_safe_harbor.slice1 import run_slice1  # noqa: E402

# This module is the ONLY place app/ is allowed to call into the core engine
# to produce a new run. It never reimplements reference validation, annotation
# parsing, interval math, or candidate rules itself (ADR-0010) -- every function
# here is a thin pass-through to the same top-level functions the CLI calls.


def trigger_baseline(
    manifest_path: Path,
    reference_key: str,
    data_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    try:
        manifest = load_reference_manifest(manifest_path)
        reference = reference_by_key(manifest, reference_key)
        paths = verify_reference_files(reference, data_dir)
        run_manifest = run_baseline(reference, paths["fasta"], paths["annotation"], output_dir)
        return {"status": "success", "run_manifest": run_manifest, "error": None}
    except (ContractError, OSError, ValueError) as exc:
        return {"status": "failed", "run_manifest": None, "error": str(exc)}


def trigger_candidate_windows(
    baseline_run_dir: Path,
    manifest_path: Path,
    reference_key: str,
    data_dir: Path,
    buffer_distance_bp: int,
    min_candidate_window_bp: int,
    output_dir: Path,
    long_interval_percentile: float = 0.95,
) -> dict[str, Any]:
    try:
        manifest = load_reference_manifest(manifest_path)
        reference = reference_by_key(manifest, reference_key)
        paths = verify_reference_files(reference, data_dir)
        run_manifest = run_slice1(
            baseline_run_dir,
            reference,
            paths["fasta"],
            buffer_distance_bp,
            min_candidate_window_bp,
            output_dir,
            long_interval_percentile,
        )
        return {"status": "success", "run_manifest": run_manifest, "error": None}
    except (ContractError, OSError, ValueError) as exc:
        return {"status": "failed", "run_manifest": None, "error": str(exc)}
