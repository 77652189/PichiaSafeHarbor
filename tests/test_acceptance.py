from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pytest

from pichia_safe_harbor.acceptance import create_acceptance_manifest
from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.pipeline import run_baseline


def _write_evidence(run_dir: Path, implementation_hash: str) -> tuple[Path, Path]:
    run_manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    statistics = json.loads((run_dir / "statistics.json").read_text(encoding="utf-8"))
    nuclear_sequences = sorted(
        seqid
        for seqid, values in statistics["sequence_statistics"].items()
        if not values["excluded_from_nuclear_statistics"]
    )
    independent = run_dir / "verification/independent_check.json"
    independent.parent.mkdir(parents=True, exist_ok=True)
    independent.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "verification_type": "independent_interval_recalculation",
                "status": "passed",
                "run_id": run_manifest["run_id"],
                "verified_intergenic_region_count": statistics["intergenic_region_count"],
                "nuclear_sequences": nuclear_sequences,
                "coverage_checks": {
                    seqid: {"conserved": True} for seqid in nuclear_sequences
                },
                "verified_artifacts": {
                    "intergenic_regions.tsv": hashlib.sha256(
                        (run_dir / "intergenic_regions.tsv").read_bytes()
                    ).hexdigest(),
                    "terminal_regions.tsv": hashlib.sha256(
                        (run_dir / "terminal_regions.tsv").read_bytes()
                    ).hexdigest(),
                },
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    tests = run_dir / "verification/test_evidence.json"
    tests.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "evidence_type": "automated_tests",
                "status": "passed",
                "command": ["python", "-m", "pytest", "-q"],
                "passed_count": 1,
                "implementation_sha256": implementation_hash,
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return independent, tests


def test_acceptance_manifest_binds_run_artifacts_and_verification(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    run_dir = tmp_path / "run"
    run = run_baseline(fixture_reference, *fixture_paths, run_dir)
    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    result = create_acceptance_manifest(run_dir, independent, tests)
    assert result["execution_status"] == "complete"
    assert result["verification_status"] == "passed"
    assert result["scientific_acceptance_status"] == "blocked"
    assert result["run_id"] == run["run_id"]
    assert set(result["artifacts"]) == set(run["artifacts"])
    assert (run_dir / "acceptance_manifest.json").is_file()


def test_acceptance_rejects_missing_or_mismatched_independent_evidence(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    run_dir = tmp_path / "run"
    run = run_baseline(fixture_reference, *fixture_paths, run_dir)
    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    independent.unlink()
    with pytest.raises(ContractError, match="missing"):
        create_acceptance_manifest(run_dir, independent, tests)
    independent.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "verification_type": "independent_interval_recalculation",
                "status": "passed",
                "run_id": "wrong",
                "verified_intergenic_region_count": 0,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ContractError, match="run_id mismatch"):
        create_acceptance_manifest(run_dir, independent, tests)


def test_acceptance_rejects_independent_artifact_identity_mismatch(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    run_dir = tmp_path / "run"
    run = run_baseline(fixture_reference, *fixture_paths, run_dir)
    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    evidence = json.loads(independent.read_text(encoding="utf-8"))
    evidence["verified_artifacts"]["intergenic_regions.tsv"] = "0" * 64
    independent.write_text(json.dumps(evidence), encoding="utf-8")
    with pytest.raises(ContractError, match="artifact identity mismatch"):
        create_acceptance_manifest(run_dir, independent, tests)


def test_acceptance_rejects_false_independent_count_or_coverage(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    run_dir = tmp_path / "run"
    run = run_baseline(fixture_reference, *fixture_paths, run_dir)
    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    evidence = json.loads(independent.read_text(encoding="utf-8"))
    evidence["verified_intergenic_region_count"] = 0
    independent.write_text(json.dumps(evidence), encoding="utf-8")
    with pytest.raises(ContractError, match="intergenic count mismatch"):
        create_acceptance_manifest(run_dir, independent, tests)

    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    evidence = json.loads(independent.read_text(encoding="utf-8"))
    evidence["coverage_checks"][evidence["nuclear_sequences"][0]]["conserved"] = False
    independent.write_text(json.dumps(evidence), encoding="utf-8")
    with pytest.raises(ContractError, match="coverage conservation mismatch"):
        create_acceptance_manifest(run_dir, independent, tests)


def test_acceptance_rejects_empty_or_malformed_test_evidence(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    run_dir = tmp_path / "run"
    run = run_baseline(fixture_reference, *fixture_paths, run_dir)
    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    evidence = json.loads(tests.read_text(encoding="utf-8"))
    evidence["passed_count"] = 0
    tests.write_text(json.dumps(evidence), encoding="utf-8")
    with pytest.raises(ContractError, match="passed_count must be positive"):
        create_acceptance_manifest(run_dir, independent, tests)

    independent, tests = _write_evidence(
        run_dir, run["software"]["implementation_sha256"]
    )
    evidence = json.loads(tests.read_text(encoding="utf-8"))
    evidence["command"] = []
    tests.write_text(json.dumps(evidence), encoding="utf-8")
    with pytest.raises(ContractError, match="command mismatch"):
        create_acceptance_manifest(run_dir, independent, tests)
