from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pichia_safe_harbor.io_utils import sha256_file  # noqa: E402

RUN_TYPE_BASELINE = "slice0_baseline"
RUN_TYPE_CANDIDATE_WINDOWS = "slice1_candidate_windows"
RUN_TYPE_UNSUPPORTED = "unsupported"

# User-facing labels for run_type -- internal identifiers above must never be
# shown to the user directly (they are project-planning codenames, not
# scientific terms).
RUN_TYPE_LABELS = {
    RUN_TYPE_BASELINE: "基因组基线数据",
    RUN_TYPE_CANDIDATE_WINDOWS: "候选安全港窗口",
    RUN_TYPE_UNSUPPORTED: "无法识别的运行",
}


def run_type_label(run_type: str) -> str:
    return RUN_TYPE_LABELS.get(run_type, run_type)

BASELINE_REQUIRED_ARTIFACTS = (
    "functional_entities.json",
    "functional_entities.tsv",
    "intergenic_regions.json",
    "intergenic_regions.tsv",
    "terminal_regions.json",
    "terminal_regions.tsv",
    "statistics.json",
    "diagnostics.json",
    "baseline_report.md",
)
CANDIDATE_WINDOWS_REQUIRED_ARTIFACTS = (
    "candidate_windows.json",
    "candidate_windows.tsv",
    "excluded_regions.json",
    "excluded_regions.tsv",
    "statistics.json",
    "candidate_windows_report.md",
)


@dataclass(frozen=True)
class RunCatalogEntry:
    """A read-only description of one run directory for Streamlit to display.

    ``displayable`` only reflects mechanical validity (recognized schema,
    hash-verified artifacts, execution_status=complete) -- it is deliberately
    independent of scientific_acceptance_status. Every page must still show
    the three status fields prominently; this project has no run today whose
    scientific_acceptance_status is anything but "blocked" (ADR-0008/0009),
    and that must never be hidden just because a run is "displayable".
    """

    run_id: str
    run_dir: Path
    run_type: str
    execution_status: str
    verification_status: str
    scientific_acceptance_status: str
    displayable: bool
    reason: str
    manifest: dict[str, Any]
    acceptance: dict[str, Any] | None


def _classify_run_type(manifest: dict[str, Any]) -> str | None:
    if "parent_run_id" in manifest and "rule_version" in manifest:
        return RUN_TYPE_CANDIDATE_WINDOWS
    if manifest.get("schema_version") == 2 and "target_strain" in manifest and "parent_run_id" not in manifest:
        return RUN_TYPE_BASELINE
    return None


def _artifacts_verified(run_dir: Path, manifest: dict[str, Any]) -> bool:
    for name, expected in manifest.get("artifacts", {}).items():
        path = run_dir / name
        if not path.is_file():
            return False
        if not isinstance(expected, dict):
            return False
        if path.stat().st_size != expected.get("size_bytes"):
            return False
        if sha256_file(path) != expected.get("sha256"):
            return False
    return True


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    return value if isinstance(value, dict) else None


def _catalog_entry_for(run_dir: Path) -> RunCatalogEntry | None:
    manifest_path = run_dir / "run_manifest.json"
    manifest = _load_json(manifest_path)
    if manifest is None:
        return None

    run_id = str(manifest.get("run_id", run_dir.name))
    execution_status = str(manifest.get("execution_status", "unknown"))
    scientific_status = str(manifest.get("scientific_acceptance_status", "unknown"))
    acceptance = _load_json(run_dir / "acceptance_manifest.json")
    verification_status = str(
        (acceptance or {}).get("verification_status", manifest.get("verification_status", "unknown"))
    )
    if acceptance is not None and acceptance.get("run_id") == run_id:
        scientific_status = str(acceptance.get("scientific_acceptance_status", scientific_status))

    run_type = _classify_run_type(manifest)
    if run_type is None:
        return RunCatalogEntry(
            run_id=run_id,
            run_dir=run_dir,
            run_type=RUN_TYPE_UNSUPPORTED,
            execution_status=execution_status,
            verification_status=verification_status,
            scientific_acceptance_status=scientific_status,
            displayable=False,
            reason="unrecognized or unsupported run_manifest schema",
            manifest=manifest,
            acceptance=acceptance,
        )

    required = BASELINE_REQUIRED_ARTIFACTS if run_type == RUN_TYPE_BASELINE else CANDIDATE_WINDOWS_REQUIRED_ARTIFACTS
    missing = [name for name in required if name not in manifest.get("artifacts", {})]
    if missing:
        displayable, reason = False, f"run_manifest is missing required artifact entries: {missing}"
    elif execution_status != "complete":
        displayable, reason = False, f"execution_status={execution_status}, not complete"
    elif not _artifacts_verified(run_dir, manifest):
        displayable, reason = False, "one or more artifacts do not match the hashes recorded in run_manifest.json"
    else:
        displayable, reason = True, "ok"

    return RunCatalogEntry(
        run_id=run_id,
        run_dir=run_dir,
        run_type=run_type,
        execution_status=execution_status,
        verification_status=verification_status,
        scientific_acceptance_status=scientific_status,
        displayable=displayable,
        reason=reason,
        manifest=manifest,
        acceptance=acceptance,
    )


def scan_run_catalog(root: Path) -> list[RunCatalogEntry]:
    """Build a read-only index of run directories directly under ``root``.

    Never modifies anything under ``root``. Directories without a
    ``run_manifest.json`` are silently skipped (they are not runs at all);
    directories with one are always included in the result, tagged
    ``displayable=False`` with a ``reason`` if they fail any check, so a
    caller can still show them on an audit/diagnostic view.
    """
    if not root.is_dir():
        return []
    entries: list[RunCatalogEntry] = []
    for child in sorted(p for p in root.iterdir() if p.is_dir()):
        entry = _catalog_entry_for(child)
        if entry is not None:
            entries.append(entry)
    return entries


def displayable_entries(entries: list[RunCatalogEntry], run_type: str | None = None) -> list[RunCatalogEntry]:
    return [
        entry
        for entry in entries
        if entry.displayable and (run_type is None or entry.run_type == run_type)
    ]
