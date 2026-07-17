from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.pipeline import _implementation_hash
from pichia_safe_harbor.transcript_probe import create_slice0b_acceptance


CHROMOSOMES = ["CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"]


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def _identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path, Path]:
    run_dir = tmp_path / "run"
    sources = tmp_path / "sources.json"
    toolchain = tmp_path / "toolchain.json"
    reference = tmp_path / "reference.json"
    acquisition = tmp_path / "acquisition.json"
    _write_json(acquisition, {"fixture": "acquisition"})
    _write_json(sources, {"fixture": True, "acquisition_evidence": acquisition.name})
    _write_json(toolchain, {"fixture": True})
    _write_json(reference, {"fixture": True})
    windows = {f"{seqid}:{window}": {"sampled_records": 1, "endpoint_affected_records": 0, "endpoint_and_multi_mapping": 0} for seqid in CHROMOSOMES for window in ("start", "middle", "end")}
    windows[f"{CHROMOSOMES[-1]}:end"]["sampled_records"] = 0
    _write_json(
        run_dir / "transcript_evidence_probe.json",
        {
            "probe_mode": "submitted_alignment_coordinate_probe",
            "records": [],
            "summary": {
                "sampled_record_count": 0,
                "sample_classification_counts": {"success": 1, "conflict": 0, "multi_mapping": 0, "unmappable": 0, "endpoint_limited": 0},
                "reported_category_counts": {"success": 1, "conflict": 0, "multi_mapping": 0, "unmappable": 0, "endpoint_limited": 0},
                "chromosomes_covered": CHROMOSOMES,
                "coordinate_compatible_chromosomes": CHROMOSOMES,
                "chromosome_coordinate_probe": "passed",
                "window_coverage": "partial",
                "coordinate_compatibility": "partial",
                "window_count": 12,
                "zero_coverage_windows": [f"{CHROMOSOMES[-1]}:end"],
                "window_statistics": windows,
                "archive_alignment_statistics": {"total_reads": 1, "primary_alignment_rows": 1, "secondary_alignment_rows": 0, "unmappable_or_not_primary_rows_estimate": 0},
                "transcription_strand_status": "unavailable; fixture",
                "controlled_unique_mapping_status": "unavailable; fixture",
                "controlled_splice_evidence_status": "unavailable; fixture",
                "true_unmappable_status": "unavailable; fixture",
                "contamination_status": "unavailable; fixture",
                "endpoint_affected_records": 0,
                "endpoint_and_multi_mapping": 0,
            }
        },
    )
    _write_json(
        run_dir / "source_applicability_matrix.json",
        {"sources": [{"eligible_as_independent_transcript_source": False, "independence_from_current_annotation": "not-passed", "license_gate": "unresolved", "strain_secondary_identifier_status": "conflict", "probe_mode": "submitted_alignment_coordinate_probe", "controlled_unique_mapping": "unavailable", "controlled_splice_evidence": "unavailable", "true_unmappable": "unavailable", "contamination": "unavailable"}]},
    )
    (run_dir / "transcript_evidence_probe.tsv").write_text("fixture\n", encoding="utf-8")
    (run_dir / "source_applicability_matrix.tsv").write_text("fixture\n", encoding="utf-8")
    (run_dir / "transcript_quality_report.md").write_text("No safe-harbor candidates or thresholds are generated\n", encoding="utf-8")
    (run_dir / "slice0b_recommendation.md").write_text("Scientific status: blocked\nDo not start full-genome reannotation or Slice 1\n", encoding="utf-8")
    names = ["transcript_evidence_probe.json", "transcript_evidence_probe.tsv", "source_applicability_matrix.json", "source_applicability_matrix.tsv", "transcript_quality_report.md", "slice0b_recommendation.md"]
    artifacts = {name: _identity(run_dir / name) for name in names}
    run_id = "slice0b-0123456789abcdef"
    _write_json(
        run_dir / "run_manifest.json",
        {
            "run_id": run_id,
            "execution_status": "complete",
            "scientific_acceptance_status": "blocked",
            "scientific_acceptance_blockers": ["fixture blocker"],
            "primary_coordinate_space": "GCA_001746955.1",
            "exact_target_strain_coordinates": False,
            "implementation_sha256": _implementation_hash(),
            "source_manifest_sha256": _identity(sources)["sha256"],
            "toolchain_sha256": _identity(toolchain)["sha256"],
            "reference_manifest_sha256": _identity(reference)["sha256"],
            "acquisition_evidence_sha256": _identity(acquisition)["sha256"],
            "artifacts": artifacts,
        },
    )
    independent = run_dir / "verification/independent_check.json"
    _write_json(
        independent,
        {
            "schema_version": 1,
            "verification_type": "slice0b_independent_summary",
            "status": "passed",
            "run_id": run_id,
            "verified_artifacts": {name: identity["sha256"] for name, identity in artifacts.items()},
            "acquisition_evidence_sha256": _identity(acquisition)["sha256"],
            "probe_mode": "submitted_alignment_coordinate_probe",
            "sampled_record_count": 0,
            "sample_classification_counts": {"success": 1, "conflict": 0, "multi_mapping": 0, "unmappable": 0, "endpoint_limited": 0},
            "chromosomes_covered": CHROMOSOMES,
            "window_count": 12,
            "zero_coverage_windows": [f"{CHROMOSOMES[-1]}:end"],
            "chromosome_coordinate_probe": "passed",
            "window_coverage": "partial",
            "coordinate_compatibility": "partial",
            "endpoint_affected_records": 0,
            "endpoint_and_multi_mapping": 0,
            "quality_gate_status": "unavailable_without_controlled_read_remapping",
            "archive_alignment_statistics": {"total_reads": 1, "primary_alignment_rows": 1, "secondary_alignment_rows": 0, "unmappable_or_not_primary_rows_estimate": 0},
            "source_independence_gate": "not-passed",
            "source_license_gate": "unresolved",
            "scientific_acceptance_status": "blocked",
            "protected_evidence": {"slice0_artifact_count": 9, "slice0a_artifact_count": 7},
        },
    )
    tests = run_dir / "verification/test_evidence.json"
    _write_json(tests, {"schema_version": 1, "evidence_type": "automated_tests", "status": "passed", "command": ["python", "-m", "pytest", "-q"], "passed_count": 52, "implementation_sha256": _implementation_hash()})
    repeatability = tmp_path / "repeatability.json"
    repeat_files = [*names, "run_manifest.json", "verification/independent_check.json", "verification/test_evidence.json"]
    peer = tmp_path / "peer"
    shutil.copytree(run_dir, peer)
    _write_json(repeatability, {"schema_version": 1, "verification_type": "slice0b_repeatability", "status": "passed", "run_id": run_id, "run_directories": [run_dir.name, peer.name], "compared_file_count": 9, "all_files_identical": True, "file_sha256": {name.replace("\\", "/"): _identity(run_dir / name)["sha256"] for name in repeat_files}})
    return run_dir, independent, tests, repeatability, peer, tmp_path, sources, toolchain, reference


def _refresh(run_dir: Path, independent: Path, repeatability: Path, peer: Path, name: str) -> None:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    run["artifacts"][name] = _identity(run_dir / name)
    _write_json(run_dir / "run_manifest.json", run)
    check = json.loads(independent.read_text(encoding="utf-8"))
    check["verified_artifacts"][name] = run["artifacts"][name]["sha256"]
    _write_json(independent, check)
    repeat = json.loads(repeatability.read_text(encoding="utf-8"))
    repeat["file_sha256"][name] = run["artifacts"][name]["sha256"]
    repeat["file_sha256"]["run_manifest.json"] = _identity(run_dir / "run_manifest.json")["sha256"]
    repeat["file_sha256"]["verification/independent_check.json"] = _identity(independent)["sha256"]
    _write_json(repeatability, repeat)
    shutil.copy2(run_dir / name, peer / name)
    shutil.copy2(run_dir / "run_manifest.json", peer / "run_manifest.json")
    shutil.copy2(independent, peer / "verification/independent_check.json")


def test_slice0b_acceptance_records_engineering_pass_and_scientific_block(tmp_path: Path) -> None:
    args = _fixture(tmp_path)
    result = create_slice0b_acceptance(*args)
    assert result["verification_status"] == "passed"
    assert result["scientific_acceptance_status"] == "blocked"
    assert result["independent_transcript_source_passed"] is False
    assert result["acquisition_evidence"]["sha256"] == _identity(tmp_path / "acquisition.json")["sha256"]


@pytest.mark.parametrize("failure", ["missing_category", "overclaimed_source", "missing_chromosome", "bad_endpoint_summary", "overclaimed_quality", "incomplete_independent", "missing_repeat_run", "missing_peer_directory", "insufficient_tests", "tampered_acquisition"])
def test_slice0b_acceptance_rejects_incomplete_or_overclaimed_evidence(tmp_path: Path, failure: str) -> None:
    run_dir, independent, tests, repeatability, peer, repo_root, sources, toolchain, reference = _fixture(tmp_path)
    if failure == "missing_category":
        path = run_dir / "transcript_evidence_probe.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        del value["summary"]["reported_category_counts"]["conflict"]
    elif failure == "overclaimed_source":
        path = run_dir / "source_applicability_matrix.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        value["sources"][0]["eligible_as_independent_transcript_source"] = True
    elif failure in {"missing_chromosome", "bad_endpoint_summary", "overclaimed_quality"}:
        path = run_dir / "transcript_evidence_probe.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        if failure == "missing_chromosome":
            value["summary"]["chromosomes_covered"] = CHROMOSOMES[:-1]
        elif failure == "bad_endpoint_summary":
            value["summary"]["endpoint_affected_records"] = 1
        else:
            value["summary"]["controlled_unique_mapping_status"] = "passed"
    elif failure == "incomplete_independent":
        path = independent
        value = json.loads(path.read_text(encoding="utf-8"))
        del value["endpoint_affected_records"]
    elif failure == "missing_repeat_run":
        path = repeatability
        value = json.loads(path.read_text(encoding="utf-8"))
        value["run_directories"] = [run_dir.name]
    elif failure == "missing_peer_directory":
        shutil.rmtree(peer)
        path = repeatability
        value = json.loads(path.read_text(encoding="utf-8"))
    elif failure == "insufficient_tests":
        path = tests
        value = json.loads(path.read_text(encoding="utf-8"))
        value["passed_count"] = 1
    else:
        path = repo_root / "acquisition.json"
        value = {"fixture": "tampered"}
    _write_json(path, value)
    if path.parent == run_dir:
        _refresh(run_dir, independent, repeatability, peer, path.name)
    with pytest.raises(ContractError):
        create_slice0b_acceptance(run_dir, independent, tests, repeatability, peer, repo_root, sources, toolchain, reference)
