from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.services.catalog import (
    RUN_TYPE_BASELINE,
    RUN_TYPE_CANDIDATE_WINDOWS,
    RUN_TYPE_UNSUPPORTED,
    scan_run_catalog,
)
from app.services.presentation import (
    candidate_view,
    load_candidate_windows,
    load_excluded_regions,
    load_statistics,
    status_banner,
)
from app.services.trigger import trigger_baseline, trigger_candidate_windows
from pichia_safe_harbor.pipeline import run_baseline
from pichia_safe_harbor.slice1 import run_slice1


def _baseline(tmp_path: Path, fixture_paths, fixture_reference, name: str = "baseline") -> Path:
    output = tmp_path / name
    run_baseline(fixture_reference, *fixture_paths, output)
    return output


def test_scan_run_catalog_classifies_baseline_and_candidate_windows_runs(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1, root / "candidates")

    entries = scan_run_catalog(root)
    by_type = {entry.run_type: entry for entry in entries}

    assert by_type[RUN_TYPE_BASELINE].displayable is True
    assert by_type[RUN_TYPE_CANDIDATE_WINDOWS].displayable is True
    assert by_type[RUN_TYPE_CANDIDATE_WINDOWS].scientific_acceptance_status == "blocked"


def test_scan_run_catalog_skips_directories_without_manifest(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    (root / "not_a_run").mkdir(parents=True)
    (root / "not_a_run" / "notes.txt").write_text("hello", encoding="utf-8")

    entries = scan_run_catalog(root)
    assert entries == []


def test_scan_run_catalog_flags_tampered_artifact_as_not_displayable(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")
    (baseline_dir / "statistics.json").write_text("{}", encoding="utf-8")

    entries = scan_run_catalog(root)
    entry = next(item for item in entries if item.run_type == RUN_TYPE_BASELINE)
    assert entry.displayable is False
    assert "hash" in entry.reason or "mismatch" in entry.reason


def test_scan_run_catalog_flags_incomplete_execution_as_not_displayable(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")
    manifest_path = baseline_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["execution_status"] = "failed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    entries = scan_run_catalog(root)
    entry = next(item for item in entries if item.run_type == RUN_TYPE_BASELINE)
    assert entry.displayable is False
    assert "execution_status" in entry.reason


def test_scan_run_catalog_marks_unrecognized_manifest_as_unsupported(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    weird = root / "weird_run"
    weird.mkdir(parents=True)
    (weird / "run_manifest.json").write_text(json.dumps({"schema_version": 99}), encoding="utf-8")

    entries = scan_run_catalog(root)
    assert len(entries) == 1
    assert entries[0].run_type == RUN_TYPE_UNSUPPORTED
    assert entries[0].displayable is False


def test_scan_run_catalog_reads_acceptance_manifest_status_when_present(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")
    manifest = json.loads((baseline_dir / "run_manifest.json").read_text(encoding="utf-8"))
    acceptance = {
        "run_id": manifest["run_id"],
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
    }
    (baseline_dir / "acceptance_manifest.json").write_text(json.dumps(acceptance), encoding="utf-8")

    entries = scan_run_catalog(root)
    entry = next(item for item in entries if item.run_type == RUN_TYPE_BASELINE)
    assert entry.verification_status == "passed"
    assert entry.acceptance == acceptance


def test_presentation_loaders_and_candidate_view(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1, root / "candidates")

    entries = scan_run_catalog(root)
    candidate_entry = next(item for item in entries if item.run_type == RUN_TYPE_CANDIDATE_WINDOWS)

    windows = load_candidate_windows(candidate_entry)
    excluded = load_excluded_regions(candidate_entry)
    statistics = load_statistics(candidate_entry)
    assert len(windows) == statistics["candidate_window_count"]
    assert len(excluded) == statistics["excluded_region_count"]

    view = candidate_view(windows[0])
    assert view["collinearity_status"] == "unavailable"
    assert "rule_flag_ids" in view

    banner = status_banner(candidate_entry)
    assert banner["scientific_acceptance_status"] == "blocked"
    assert banner["exact_target_strain_coordinates"] is False


def test_presentation_loaders_reject_wrong_run_type(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    root = tmp_path / "runs"
    root.mkdir()
    _baseline(tmp_path, fixture_paths, fixture_reference, name="runs/baseline")

    entries = scan_run_catalog(root)
    baseline_entry = next(item for item in entries if item.run_type == RUN_TYPE_BASELINE)

    with pytest.raises(ValueError, match="not a candidate-windows run"):
        load_candidate_windows(baseline_entry)


def test_trigger_baseline_reports_failure_without_raising(tmp_path: Path) -> None:
    result = trigger_baseline(
        tmp_path / "missing_manifest.json", "strain-b", tmp_path / "data", tmp_path / "out"
    )
    assert result["status"] == "failed"
    assert result["run_manifest"] is None
    assert result["error"]


def test_trigger_candidate_windows_reports_failure_without_raising(tmp_path: Path) -> None:
    result = trigger_candidate_windows(
        tmp_path / "missing_baseline",
        tmp_path / "missing_manifest.json",
        "strain-b",
        tmp_path / "data",
        0,
        1,
        tmp_path / "out",
    )
    assert result["status"] == "failed"
    assert result["run_manifest"] is None
    assert result["error"]
