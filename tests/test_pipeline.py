from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.pipeline import run_baseline


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_for_paths(reference: dict, fasta_path: Path, gff_path: Path) -> dict:
    updated = deepcopy(reference)
    updated["files"] = {
        "fasta": {"sha256": _hash(fasta_path), "size_bytes": fasta_path.stat().st_size},
        "annotation": {"sha256": _hash(gff_path), "size_bytes": gff_path.stat().st_size},
    }
    return updated


def test_baseline_is_complete_and_deterministic(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first_manifest = run_baseline(fixture_reference, *fixture_paths, first)
    second_manifest = run_baseline(fixture_reference, *fixture_paths, second)
    assert first_manifest["execution_status"] == "complete"
    assert first_manifest["verification_status"] == "not_run"
    assert first_manifest["scientific_acceptance_status"] == "blocked"
    assert "status" not in first_manifest
    assert first_manifest["run_id"] == second_manifest["run_id"]
    for name in (
        "functional_entities.json",
        "functional_entities.tsv",
        "intergenic_regions.json",
        "intergenic_regions.tsv",
        "terminal_regions.json",
        "terminal_regions.tsv",
        "statistics.json",
        "diagnostics.json",
        "baseline_report.md",
        "run_manifest.json",
    ):
        assert _hash(first / name) == _hash(second / name)
    manifest = json.loads((first / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["target_strain"] == "Strain-T"
    assert manifest["primary_assembly"] == "GCA_001746955.1"
    assert manifest["exact_target_strain_coordinates"] is False
    assert "no candidate windows" in manifest["analysis_scope"]
    assert len(manifest["software"]["implementation_sha256"]) == 64
    entities = json.loads((first / "functional_entities.json").read_text(encoding="utf-8"))
    assert len(entities) == 11
    assert next(item for item in entities if item["entity_id"] == "g2")["partial"] is True
    report = (first / "baseline_report.md").read_text(encoding="utf-8")
    assert "## Annotation applicability" in report
    assert "threshold inference: unavailable" in report
    assert "Scientific acceptance: blocked" in report


def test_contract_failure_publishes_no_partial_run(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    fasta_path, gff_path = fixture_paths
    bad_gff = tmp_path / "bad.gff3"
    bad_gff.write_text(
        gff_path.read_text(encoding="utf-8")
        + "chr1\tfixture\tgene\t99\t101\t.\t+\t.\tID=out_of_bounds\n",
        encoding="utf-8",
    )
    output = tmp_path / "failed"
    with pytest.raises(ContractError, match="out of bounds"):
        run_baseline(_reference_for_paths(fixture_reference, fasta_path, bad_gff), fasta_path, bad_gff, output)
    assert not output.exists()


def test_conflicting_duplicate_id_fails(fixture_paths, fixture_reference, tmp_path: Path) -> None:
    fasta_path, gff_path = fixture_paths
    bad_gff = tmp_path / "duplicate.gff3"
    bad_gff.write_text(
        gff_path.read_text(encoding="utf-8")
        + "chr1\tfixture\tgene\t91\t95\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    with pytest.raises(ContractError, match="conflicting duplicate"):
        run_baseline(
            _reference_for_paths(fixture_reference, fasta_path, bad_gff),
            fasta_path,
            bad_gff,
            tmp_path / "failed",
        )


def test_segmented_child_id_is_allowed_and_merged(
    fixture_paths, fixture_reference, tmp_path: Path
) -> None:
    fasta_path, gff_path = fixture_paths
    segmented = tmp_path / "segmented.gff3"
    segmented.write_text(
        gff_path.read_text(encoding="utf-8")
        + "chr2\tfixture\tCDS\t57\t58\t.\t+\t0\tID=split_cds;Parent=missing_split_parent\n"
        + "chr2\tfixture\tCDS\t60\t62\t.\t+\t2\tID=split_cds;Parent=missing_split_parent\n",
        encoding="utf-8",
    )
    output = tmp_path / "segmented"
    run_baseline(
        _reference_for_paths(fixture_reference, fasta_path, segmented),
        fasta_path,
        segmented,
        output,
    )
    diagnostics = json.loads((output / "diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["duplicate_ids"][-1]["id"] == "split_cds"


def test_core_entry_rejects_actual_hash_mismatch(
    fixture_paths, fixture_reference, tmp_path: Path
) -> None:
    fasta_path, gff_path = fixture_paths
    changed = tmp_path / "changed.fna"
    changed.write_bytes(fasta_path.read_bytes().replace(b"A", b"C", 1))
    output = tmp_path / "failed"
    with pytest.raises(ContractError, match="sha256 mismatch"):
        run_baseline(fixture_reference, changed, gff_path, output)
    assert not output.exists()


def test_run_id_changes_with_actual_input_and_sequence_classification(
    fixture_paths, fixture_reference, tmp_path: Path
) -> None:
    fasta_path, gff_path = fixture_paths
    baseline = run_baseline(fixture_reference, fasta_path, gff_path, tmp_path / "baseline")

    changed_fasta = tmp_path / "changed.fna"
    changed_fasta.write_bytes(fasta_path.read_bytes().replace(b"A", b"C", 1))
    changed_reference = _reference_for_paths(fixture_reference, changed_fasta, gff_path)
    changed = run_baseline(changed_reference, changed_fasta, gff_path, tmp_path / "changed")
    assert changed["run_id"] != baseline["run_id"]

    classification_reference = deepcopy(fixture_reference)
    classification_reference["sequence_classes"]["scaffold"] = "mitochondrial"
    classified = run_baseline(
        classification_reference, fasta_path, gff_path, tmp_path / "classified"
    )
    assert classified["run_id"] != baseline["run_id"]
