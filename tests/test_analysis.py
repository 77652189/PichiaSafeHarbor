from __future__ import annotations

from collections import defaultdict

import pytest

from pichia_safe_harbor.analysis import (
    build_inventory,
    merge_boundaries,
    map_annotation_sequences,
    normalize_entities,
    validate_coordinates,
)
from pichia_safe_harbor.models import FunctionalEntity
from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.parsers import read_fasta_index, read_gff3


def _analyze(fixture_paths, fixture_reference):
    fasta_path, gff_path = fixture_paths
    fasta = read_fasta_index(fasta_path)
    features, _, _ = read_gff3(gff_path)
    features, _ = map_annotation_sequences(features, fasta, {})
    validate_coordinates(features, fasta)
    entities, diagnostics = normalize_entities(features)
    inventory, terminal, stats = build_inventory(
        fasta,
        entities,
        fixture_reference["sequence_classes"],
        fixture_reference["sequence_endpoints"],
    )
    return fasta, entities, diagnostics, inventory, terminal, stats


def test_inventory_covers_orientations_parent_graph_and_exclusions(
    fixture_paths, fixture_reference
) -> None:
    _, entities, diagnostics, inventory, terminal, stats = _analyze(
        fixture_paths, fixture_reference
    )
    chr1 = [item for item in inventory if item.seqid == "chr1"]
    assert [(item.start, item.end, item.orientation) for item in chr1] == [
        (10, 20, "convergent"),
        (30, 40, "divergent"),
        (50, 60, "tandem"),
        (70, 80, "unknown"),
    ]
    assert not {"mito", "scaffold"} & {item.seqid for item in inventory}
    assert stats["sequence_statistics"]["mito"]["excluded_from_nuclear_statistics"]
    assert stats["sequence_statistics"]["scaffold"]["excluded_from_nuclear_statistics"]
    assert stats["sequence_statistics"]["chr_empty"]["intergenic_region_count"] == 0
    assert set(stats["orientation_counts"]) == {
        "convergent",
        "tandem",
        "divergent",
        "unknown",
    }
    assert any(item.side == "whole_sequence_no_features" for item in terminal)
    chr1_3p = next(item for item in terminal if item.seqid == "chr1" and item.side == "3p")
    assert chr1_3p.adjacent_end_completeness == "partial"
    assert chr1_3p.interpretable_as_true_telomere_distance is False
    assert any(item.entity_id == "orphan_cds" for item in entities)
    assert diagnostics.missing_parents == [
        {"line": 22, "feature_id": "orphan_cds", "parent": "missing_transcript"}
    ]
    assert diagnostics.nested


def test_union_has_no_intron_or_nested_false_gaps(fixture_paths, fixture_reference) -> None:
    _, _, _, inventory, _, _ = _analyze(fixture_paths, fixture_reference)
    chr2 = [(item.start, item.end) for item in inventory if item.seqid == "chr2"]
    assert chr2 == [(30, 34), (40, 49)]


def test_statistics_partition_each_nuclear_sequence(fixture_paths, fixture_reference) -> None:
    fasta, entities, _, inventory, terminal, _ = _analyze(fixture_paths, fixture_reference)
    by_entity = defaultdict(list)
    for entity in entities:
        by_entity[entity.seqid].append(entity)
    for seqid in ("chr1", "chr2"):
        ordered = sorted(by_entity[seqid], key=lambda item: (item.start, item.end))
        merged = []
        for item in ordered:
            if merged and item.start < merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], item.end)
            else:
                merged.append([item.start, item.end])
        covered = sum(end - start for start, end in merged)
        gaps = sum(item.length for item in inventory if item.seqid == seqid)
        ends = sum(item.length for item in terminal if item.seqid == seqid)
        assert covered + gaps + ends == fasta[seqid].length


def test_unknown_sequence_classification_fails(fixture_paths, fixture_reference) -> None:
    fasta_path, gff_path = fixture_paths
    fasta = read_fasta_index(fasta_path)
    features, _, _ = read_gff3(gff_path)
    entities, _ = normalize_entities(features)
    classes = dict(fixture_reference["sequence_classes"])
    classes.pop("scaffold")
    with pytest.raises(ContractError, match="unclassified"):
        build_inventory(fasta, entities, classes, fixture_reference["sequence_endpoints"])


def test_merged_cluster_uses_gap_adjacent_boundary_strand() -> None:
    cluster = merge_boundaries(
        [
            FunctionalEntity("outer", "chr", 0, 20, "+", "gene"),
            FunctionalEntity("nested", "chr", 5, 10, "-", "gene"),
        ]
    )[0]
    assert cluster.right_boundary_entity_ids == ("outer",)
    assert cluster.right_boundary_strand == "+"


def test_partial_and_range_attributes_propagate_to_intergenic_boundaries(
    fixture_paths, fixture_reference
) -> None:
    _, entities, _, inventory, _, stats = _analyze(fixture_paths, fixture_reference)
    g1 = next(item for item in entities if item.entity_id == "g1")
    g2 = next(item for item in entities if item.entity_id == "g2")
    assert g1.partial is False
    assert g2.partial is True
    assert g2.start_range == (".", "21")
    assert g2.end_range == ("30", ".")
    assert g2.five_prime_confidence == "uncertain"
    assert g2.three_prime_confidence == "uncertain"
    first = next(item for item in inventory if item.seqid == "chr1")
    assert first.left_boundary_confidence == "high"
    assert first.right_boundary_confidence == "uncertain"
    assert stats["high_confidence_intergenic"]["count"] <= stats["intergenic_region_count"]
