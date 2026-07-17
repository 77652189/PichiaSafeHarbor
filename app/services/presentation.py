from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from .catalog import RUN_TYPE_BASELINE, RUN_TYPE_CANDIDATE_WINDOWS, RunCatalogEntry


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_statistics(entry: RunCatalogEntry) -> dict[str, Any]:
    return _read_json(entry.run_dir / "statistics.json")


def load_candidate_windows(entry: RunCatalogEntry) -> list[dict[str, Any]]:
    if entry.run_type != RUN_TYPE_CANDIDATE_WINDOWS:
        raise ValueError(f"{entry.run_id} is not a candidate-windows run")
    return _read_json(entry.run_dir / "candidate_windows.json")


def load_excluded_regions(entry: RunCatalogEntry) -> list[dict[str, Any]]:
    if entry.run_type != RUN_TYPE_CANDIDATE_WINDOWS:
        raise ValueError(f"{entry.run_id} is not a candidate-windows run")
    return _read_json(entry.run_dir / "excluded_regions.json")


def load_intergenic_regions(entry: RunCatalogEntry) -> list[dict[str, Any]]:
    if entry.run_type != RUN_TYPE_BASELINE:
        raise ValueError(f"{entry.run_id} is not a baseline run")
    return _read_json(entry.run_dir / "intergenic_regions.json")


def load_terminal_regions(entry: RunCatalogEntry) -> list[dict[str, Any]]:
    if entry.run_type != RUN_TYPE_BASELINE:
        raise ValueError(f"{entry.run_id} is not a baseline run")
    return _read_json(entry.run_dir / "terminal_regions.json")


def load_report_text(entry: RunCatalogEntry) -> str:
    name = "baseline_report.md" if entry.run_type == RUN_TYPE_BASELINE else "candidate_windows_report.md"
    path = entry.run_dir / name
    return path.read_text(encoding="utf-8") if path.is_file() else ""


def candidate_view(record: dict[str, Any]) -> dict[str, Any]:
    """Read-only display model for one candidate window.

    Formats existing fields only; ``collinearity_status`` is always
    ``unavailable`` because Slice 2 (Strain-B-Strain-C collinearity) is not
    implemented (ADR-0011) -- this must never be fabricated or silently
    dropped from the view.
    """
    view = dict(record)
    view["collinearity_status"] = "unavailable"
    view["rule_flag_ids"] = ", ".join(flag["rule_id"] for flag in record.get("rule_flags", []))
    return view


def status_banner(entry: RunCatalogEntry) -> dict[str, Any]:
    """The status facts every page must show prominently for a given run."""
    manifest = entry.manifest
    return {
        "run_id": entry.run_id,
        "run_type": entry.run_type,
        "execution_status": entry.execution_status,
        "verification_status": entry.verification_status,
        "scientific_acceptance_status": entry.scientific_acceptance_status,
        "target_strain": manifest.get("target_strain", "Strain-T"),
        "primary_reference_strain": manifest.get("primary_reference_strain", "Strain-B"),
        "exact_target_strain_coordinates": manifest.get("exact_target_strain_coordinates", False),
        "run_purpose": manifest.get("run_purpose"),
        "scientific_acceptance_blockers": manifest.get("scientific_acceptance_blockers", []),
    }
