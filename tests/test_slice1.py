from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from pichia_safe_harbor.candidates import classify_candidate_windows
from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.models import ExclusionZone, IntergenicRegion
from pichia_safe_harbor.pipeline import run_baseline
from pichia_safe_harbor.slice1 import create_slice1_acceptance, run_slice1


def _region(seqid: str = "chr1", start: int = 0, end: int = 100, orientation: str = "convergent") -> IntergenicRegion:
    return IntergenicRegion(
        region_id=f"{seqid}:intergenic:000001",
        seqid=seqid,
        start=start,
        end=end,
        length=end - start,
        left_entity_ids=("gL",),
        right_entity_ids=("gR",),
        left_strand="+",
        right_strand="-",
        orientation=orientation,
        gap_length_bp=end - start,
        left_boundary_confidence="high",
        right_boundary_confidence="high",
    )


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _baseline(tmp_path: Path, fixture_paths, fixture_reference, name: str = "baseline") -> Path:
    output = tmp_path / name
    run_baseline(fixture_reference, *fixture_paths, output)
    return output


def _by_region(records: list[dict], seqid: str, start: int, end: int) -> dict:
    return next(r for r in records if r["seqid"] == seqid and r["start"] == start and r["end"] == end)


def _independent_check(slice1_dir: Path) -> dict:
    manifest = json.loads((slice1_dir / "run_manifest.json").read_text(encoding="utf-8"))
    statistics = json.loads((slice1_dir / "statistics.json").read_text(encoding="utf-8"))
    return {
        "schema_version": 1,
        "verification_type": "slice1_independent_recalculation",
        "status": "passed",
        "run_id": manifest["run_id"],
        "verified_artifacts": {name: identity["sha256"] for name, identity in manifest["artifacts"].items()},
        "verified_candidate_window_count": statistics["candidate_window_count"],
        "verified_excluded_region_count": statistics["excluded_region_count"],
        "verified_excluded_by_rule": statistics["excluded_by_rule"],
    }


def _test_evidence(slice1_dir: Path) -> dict:
    manifest = json.loads((slice1_dir / "run_manifest.json").read_text(encoding="utf-8"))
    return {
        "schema_version": 1,
        "evidence_type": "automated_tests",
        "status": "passed",
        "command": ["python", "-m", "pytest", "-q"],
        "passed_count": 104,
        "implementation_sha256": manifest["software"]["implementation_sha256"],
    }


def _write_json(path: Path, value: dict) -> Path:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")
    return path


def test_boundary_confidence_gate_excludes_regions_next_to_partial_gene(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    manifest = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], buffer_distance_bp=0, min_candidate_window_bp=1, output_dir=output)

    assert manifest["execution_status"] == "complete"
    assert manifest["verification_status"] == "not_run"
    assert manifest["scientific_acceptance_status"] == "blocked"
    assert manifest["parent_run_id"] == json.loads((baseline_dir / "run_manifest.json").read_text(encoding="utf-8"))["run_id"]

    excluded = json.loads((output / "excluded_regions.json").read_text(encoding="utf-8"))
    windows = json.loads((output / "candidate_windows.json").read_text(encoding="utf-8"))

    # chr1 g1|g2 and g2|g3 are adjacent to the partial gene g2 (start_range and end_range set)
    # and must be excluded even though buffer/min-length would otherwise allow them.
    g1_g2 = _by_region(excluded, "chr1", 10, 20)
    assert g1_g2["rule_id"] == "boundary_confidence_insufficient"
    g2_g3 = _by_region(excluded, "chr1", 30, 40)
    assert g2_g3["rule_id"] == "boundary_confidence_insufficient"

    # chr1 g3|g4 (tandem, both high confidence) and g4|g5 (unknown orientation, both high confidence)
    # survive the boundary gate and become flagged candidate windows.
    g3_g4 = _by_region(windows, "chr1", 50, 60)
    assert g3_g4["structural_tier"] == "flagged"
    assert g3_g4["rule_flags"][0]["rule_id"] == "neighbor_orientation_tandem"
    g4_g5 = _by_region(windows, "chr1", 70, 80)
    assert g4_g5["structural_tier"] == "flagged"
    assert g4_g5["rule_flags"][0]["rule_id"] == "neighbor_orientation_unknown"

    total_input = json.loads((baseline_dir / "intergenic_regions.json").read_text(encoding="utf-8"))
    assert len(windows) + len(excluded) == len(total_input)


def test_minimum_window_length_excludes_short_remainder(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], buffer_distance_bp=0, min_candidate_window_bp=11, output_dir=output)

    excluded = json.loads((output / "excluded_regions.json").read_text(encoding="utf-8"))
    g3_g4 = _by_region(excluded, "chr1", 50, 60)
    assert g3_g4["rule_id"] == "below_minimum_window_length"


def test_buffer_distance_shrinks_candidate_window(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], buffer_distance_bp=3, min_candidate_window_bp=1, output_dir=output)

    windows = json.loads((output / "candidate_windows.json").read_text(encoding="utf-8"))
    window = _by_region(windows, "chr1", 53, 57)
    assert window["length"] == 4


def test_deterministic_across_repeated_runs(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    first = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "first")
    second = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "second")
    assert first["run_id"] == second["run_id"]
    for name in (
        "candidate_windows.json",
        "candidate_windows.tsv",
        "excluded_regions.json",
        "excluded_regions.tsv",
        "statistics.json",
        "candidate_windows_report.md",
        "run_manifest.json",
    ):
        assert _hash(tmp_path / "first" / name) == _hash(tmp_path / "second" / name)


def test_run_id_changes_with_rule_params(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    first = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "first")
    second = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 5, 1, tmp_path / "second")
    assert first["run_id"] != second["run_id"]


def test_run_manifest_declares_placeholder_and_scientific_blockers(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    manifest = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "slice1")
    assert manifest["run_purpose"] == "engine_readiness_test_not_scientific_output"
    assert manifest["target_strain"] == "Strain-T"
    assert manifest["exact_target_strain_coordinates"] is False
    assert "rule_parameters_are_illustrative_placeholders_pending_threshold_adr" in manifest["scientific_acceptance_blockers"]
    assert "no_target_strain_native_sequencing_data_used" in manifest["scientific_acceptance_blockers"]
    report = (tmp_path / "slice1" / "candidate_windows_report.md").read_text(encoding="utf-8")
    assert "NOT SCIENTIFIC OUTPUT" in report
    assert "ILLUSTRATIVE PLACEHOLDERS" in report


def test_rejects_baseline_that_is_not_complete(tmp_path: Path, fixture_paths, fixture_reference) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    manifest_path = baseline_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["execution_status"] = "in_progress"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ContractError, match="execution_status is not complete"):
        run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "slice1")


def test_rejects_missing_intergenic_regions_file(tmp_path: Path, fixture_paths, fixture_reference) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    (baseline_dir / "intergenic_regions.json").unlink()
    with pytest.raises(ContractError, match="missing intergenic_regions.json"):
        run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,tmp_path / "slice1")


def test_rejects_existing_output_dir(tmp_path: Path, fixture_paths, fixture_reference) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    output.mkdir()
    with pytest.raises(ContractError, match="output directory already exists"):
        run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,output)


def test_candidate_window_includes_actual_sequence_from_fasta(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1, output)

    windows = json.loads((output / "candidate_windows.json").read_text(encoding="utf-8"))
    window = _by_region(windows, "chr1", 50, 60)
    # the chr1 fixture sequence is a uniform run of "A", so this is a strong check
    # that extraction used the right seqid/coordinates rather than the wrong slice
    assert window["sequence"] == "A" * window["length"]

    fasta_text = (output / "candidate_windows.fasta").read_text(encoding="utf-8")
    assert f">{window['candidate_id']}" in fasta_text
    assert f"start={window['start']}" in fasta_text
    assert "A" * window["length"] in fasta_text


def test_run_slice1_rejects_fasta_not_matching_baseline(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    other_fasta = tmp_path / "different.fna"
    other_fasta.write_text(">chr1 different sequence\nCCCCCCCCCC\n", encoding="utf-8")

    with pytest.raises(ContractError, match="does not match the one used"):
        run_slice1(baseline_dir, fixture_reference, other_fasta, 0, 1, tmp_path / "slice1")


def test_candidate_windows_fasta_artifact_is_recorded_in_manifest(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    output = tmp_path / "slice1"
    manifest = run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1, output)

    assert "candidate_windows.fasta" in manifest["artifacts"]
    recorded = manifest["artifacts"]["candidate_windows.fasta"]
    actual_path = output / "candidate_windows.fasta"
    assert actual_path.is_file()
    assert recorded["size_bytes"] == actual_path.stat().st_size
    assert recorded["sha256"] == _hash(actual_path)


def test_classify_candidate_windows_rejects_invalid_rule_params() -> None:
    region = IntergenicRegion(
        region_id="chr1:intergenic:000001",
        seqid="chr1",
        start=0,
        end=10,
        length=10,
        left_entity_ids=(),
        right_entity_ids=(),
        left_strand="+",
        right_strand="-",
        orientation="convergent",
        gap_length_bp=10,
        left_boundary_confidence="high",
        right_boundary_confidence="high",
    )
    with pytest.raises(ValueError, match="buffer_distance_bp"):
        classify_candidate_windows([region], -1, 1)
    with pytest.raises(ValueError, match="min_candidate_window_bp"):
        classify_candidate_windows([region], 0, 0)


def test_convergent_high_confidence_region_is_clean_with_no_flags() -> None:
    region = IntergenicRegion(
        region_id="chr1:intergenic:000001",
        seqid="chr1",
        start=0,
        end=10,
        length=10,
        left_entity_ids=("g1",),
        right_entity_ids=("g2",),
        left_strand="+",
        right_strand="-",
        orientation="convergent",
        gap_length_bp=10,
        left_boundary_confidence="high",
        right_boundary_confidence="high",
    )
    windows, excluded, summary = classify_candidate_windows([region], 0, 1, long_interval_percentile=0.95)
    assert not excluded
    assert windows[0].structural_tier == "convergent_clean"
    assert windows[0].rule_flags == ()
    assert summary["candidate_window_count"] == 1


def test_create_slice1_acceptance_passes_with_valid_evidence(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    independent_path = _write_json(tmp_path / "independent.json", _independent_check(slice1_dir))
    test_evidence_path = _write_json(tmp_path / "tests.json", _test_evidence(slice1_dir))

    result = create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)

    statistics = json.loads((slice1_dir / "statistics.json").read_text(encoding="utf-8"))
    assert result["execution_status"] == "complete"
    assert result["verification_status"] == "passed"
    assert result["scientific_acceptance_status"] == "blocked"
    assert result["run_purpose"] == "engine_readiness_test_not_scientific_output"
    assert result["candidate_summary"]["candidate_window_count"] == statistics["candidate_window_count"]
    assert result["candidate_summary"]["candidate_window_count"] > 0
    manifest_path = slice1_dir / "acceptance_manifest.json"
    assert manifest_path.is_file()
    assert json.loads(manifest_path.read_text(encoding="utf-8")) == result


def test_create_slice1_acceptance_rejects_run_id_mismatch(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    independent = _independent_check(slice1_dir)
    independent["run_id"] = "slice1-0000000000000000"
    independent_path = _write_json(tmp_path / "independent.json", independent)
    test_evidence_path = _write_json(tmp_path / "tests.json", _test_evidence(slice1_dir))

    with pytest.raises(ContractError, match="run_id mismatch"):
        create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)


def test_create_slice1_acceptance_rejects_candidate_count_mismatch(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    independent = _independent_check(slice1_dir)
    independent["verified_candidate_window_count"] += 1
    independent_path = _write_json(tmp_path / "independent.json", independent)
    test_evidence_path = _write_json(tmp_path / "tests.json", _test_evidence(slice1_dir))

    with pytest.raises(ContractError, match="candidate window count mismatch"):
        create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)


def test_create_slice1_acceptance_rejects_tampered_artifact(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    independent_path = _write_json(tmp_path / "independent.json", _independent_check(slice1_dir))
    test_evidence_path = _write_json(tmp_path / "tests.json", _test_evidence(slice1_dir))
    (slice1_dir / "candidate_windows.json").write_text("[]", encoding="utf-8")

    with pytest.raises(ContractError, match="artifact identity mismatch"):
        create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)


def test_create_slice1_acceptance_rejects_failed_test_evidence(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    independent_path = _write_json(tmp_path / "independent.json", _independent_check(slice1_dir))
    test_evidence = _test_evidence(slice1_dir)
    test_evidence["status"] = "failed"
    test_evidence_path = _write_json(tmp_path / "tests.json", test_evidence)

    with pytest.raises(ContractError, match="automated tests did not pass"):
        create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)


def test_create_slice1_acceptance_rejects_incomplete_run(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    baseline_dir = _baseline(tmp_path, fixture_paths, fixture_reference)
    slice1_dir = tmp_path / "slice1"
    run_slice1(baseline_dir, fixture_reference, fixture_paths[0], 0, 1,slice1_dir)

    manifest_path = slice1_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["execution_status"] = "in_progress"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    independent_path = _write_json(tmp_path / "independent.json", _independent_check(slice1_dir))
    test_evidence_path = _write_json(tmp_path / "tests.json", _test_evidence(slice1_dir))

    with pytest.raises(ContractError, match="execution_status=complete"):
        create_slice1_acceptance(slice1_dir, independent_path, test_evidence_path)


def test_no_zones_produces_single_unsplit_window() -> None:
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 5, extra_exclusion_zones=[])
    assert not excluded
    assert len(windows) == 1
    assert (windows[0].start, windows[0].end) == (0, 100)
    assert windows[0].split_index == 1
    assert windows[0].split_count == 1
    assert windows[0].left_entity_ids == ("gL",)
    assert windows[0].right_entity_ids == ("gR",)


def test_single_interior_zone_splits_region_into_two_windows() -> None:
    zone = ExclusionZone(seqid="chr1", start=40, end=60, rule_id="repeat_element", reason="repeat track overlap")
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 5, extra_exclusion_zones=[zone])

    assert summary["candidate_window_count"] == 2
    ordered = sorted(windows, key=lambda item: item.start)
    first, second = ordered
    assert (first.start, first.end) == (0, 40)
    assert (second.start, second.end) == (60, 100)
    assert first.split_index == 1 and first.split_count == 2
    assert second.split_index == 2 and second.split_count == 2
    # only the window touching the original region boundary keeps that side's neighbor gene
    assert first.left_entity_ids == ("gL",) and first.right_entity_ids == ()
    assert second.left_entity_ids == () and second.right_entity_ids == ("gR",)

    zone_cut = next(item for item in excluded if item.rule_id == "repeat_element")
    assert (zone_cut.start, zone_cut.end) == (40, 60)
    assert zone_cut.parent_region_id == "chr1:intergenic:000001"


def test_two_interior_zones_split_region_into_three_windows() -> None:
    zones = [
        ExclusionZone(seqid="chr1", start=20, end=30, rule_id="repeat_element", reason="repeat a"),
        ExclusionZone(seqid="chr1", start=70, end=80, rule_id="repeat_element", reason="repeat b"),
    ]
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 5, extra_exclusion_zones=zones)

    assert summary["candidate_window_count"] == 3
    spans = sorted((item.start, item.end) for item in windows)
    assert spans == [(0, 20), (30, 70), (80, 100)]
    assert {item.rule_id for item in excluded} == {"repeat_element"}
    assert len([item for item in excluded if item.rule_id == "repeat_element"]) == 2


def test_zone_leaves_only_short_fragments_which_are_excluded() -> None:
    zone = ExclusionZone(seqid="chr1", start=5, end=95, rule_id="repeat_element", reason="repeat covers middle")
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 10, extra_exclusion_zones=[zone])

    assert windows == []
    reasons_by_rule = _count_by_rule(excluded)
    assert reasons_by_rule["repeat_element"] == 1
    assert reasons_by_rule["below_minimum_window_length"] == 2
    for item in excluded:
        if item.rule_id == "below_minimum_window_length":
            assert item.length < 10


def test_zone_fully_covering_span_leaves_no_fragments() -> None:
    zone = ExclusionZone(seqid="chr1", start=0, end=100, rule_id="repeat_element", reason="repeat covers everything")
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 5, extra_exclusion_zones=[zone])

    assert windows == []
    assert len(excluded) == 1
    assert excluded[0].rule_id == "repeat_element"


def test_zone_outside_candidate_span_has_no_effect() -> None:
    zone = ExclusionZone(seqid="chr1", start=200, end=210, rule_id="repeat_element", reason="unrelated region")
    windows, excluded, summary = classify_candidate_windows([_region()], 0, 5, extra_exclusion_zones=[zone])

    assert not excluded
    assert len(windows) == 1
    assert windows[0].split_count == 1


def _count_by_rule(excluded: list) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in excluded:
        counts[item.rule_id] = counts.get(item.rule_id, 0) + 1
    return counts
