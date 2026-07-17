from __future__ import annotations

import random

from pichia_safe_harbor.mapping_probe import map_feature_boundaries, reverse_complement, select_probe_windows
from pichia_safe_harbor.models import Feature


def _sequence(seed: int, length: int = 500) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _feature(seqid: str = "old1", start: int = 150, end: int = 260, strand: str = "+", identifier: str = "gene1") -> Feature:
    return Feature(seqid, "fixture", "gene", start, end, strand, {"ID": (identifier,), "Name": (identifier,), "locus_tag": (identifier,)}, 1)


def test_mapping_probe_classifies_consistent_offset_and_reverse() -> None:
    source = _sequence(1)
    feature = _feature()
    target = _sequence(2, 80) + source + _sequence(3, 80)
    mapped = map_feature_boundaries(feature, source, {"main": target})
    assert mapped["mapping_status"] == "consistent"
    assert mapped["assembly_orientation"] == "+"
    assert mapped["target_start"] == feature.start + 80
    reverse_target = _sequence(4, 60) + reverse_complement(source) + _sequence(5, 60)
    reverse = map_feature_boundaries(feature, source, {"main": reverse_target})
    assert reverse["mapping_status"] == "consistent"
    assert reverse["assembly_orientation"] == "-"
    assert reverse["mapped_strand"] == "-"


def test_mapping_probe_classifies_missing_duplicate_and_conflict() -> None:
    source = _sequence(10)
    feature = _feature()
    missing = map_feature_boundaries(feature, source, {"main": _sequence(11)})
    assert missing["mapping_status"] == "unmappable"
    duplicated = map_feature_boundaries(feature, source, {"main": source + _sequence(12, 20) + source})
    assert duplicated["mapping_status"] == "uncertain"
    flank = 50
    start_anchor = source[feature.start - flank : feature.start + flank + 1]
    end_center = feature.end - 1
    end_anchor = source[end_center - flank : end_center + flank + 1]
    conflict = map_feature_boundaries(feature, source, {"main1": _sequence(13, 80) + start_anchor + _sequence(14, 80), "main2": _sequence(15, 80) + end_anchor + _sequence(16, 80)})
    assert conflict["mapping_status"] == "conflict"


def test_probe_selection_covers_all_four_chromosomes() -> None:
    features = []
    for chromosome in range(1, 5):
        for index in range(10):
            features.append(_feature(f"old{chromosome}", 100 + index * 300, 200 + index * 300, identifier=f"g{chromosome}_{index}"))
    selected = select_probe_windows(features)
    assert {feature.seqid for _, _, feature in selected} == {"old1", "old2", "old3", "old4"}
    assert {label for _, label, _ in selected} == {"first", "middle", "last"}


def test_mapping_probe_does_not_overclaim_gene_adjacency(tmp_path) -> None:
    from pichia_safe_harbor.mapping_probe import run_mapping_probe

    source = _sequence(40, 2400)
    source_fasta = tmp_path / "source.fna"
    target_fasta = tmp_path / "target.fna"
    source_gff = tmp_path / "source.gff3"
    source_fasta.write_text(f">old1\n{source}\n", encoding="utf-8")
    target_fasta.write_text(f">main\n{_sequence(41, 100)}{source}\n", encoding="utf-8")
    genes = []
    for index in range(9):
        start = 151 + index * 220
        end = start + 120
        genes.append(f"old1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=g{index};Name=g{index};locus_tag=g{index}")
    source_gff.write_text("##gff-version 3\n" + "\n".join(genes) + "\n", encoding="utf-8")
    records, summary = run_mapping_probe(source_fasta, source_gff, target_fasta, {"main"})
    assert summary["order_status_counts"] == {"preserved": 3}
    assert summary["adjacency_status_counts"] == {"unavailable": 3}
    assert {record["adjacency_status"] for record in records} == {"unavailable"}
