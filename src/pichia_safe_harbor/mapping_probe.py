from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from .models import Feature
from .parsers import read_fasta_sequences, read_gff3


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _hits(anchor: str, targets: dict[str, str], max_hits: int = 3) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    reverse = reverse_complement(anchor)
    for seqid, sequence in sorted(targets.items()):
        for query, orientation in ((anchor, "+"), (reverse, "-")):
            start = 0
            while len(hits) < max_hits:
                position = sequence.find(query, start)
                if position < 0:
                    break
                item = {"seqid": seqid, "anchor_start": position, "orientation": orientation}
                if item not in hits:
                    hits.append(item)
                start = position + 1
        if len(hits) >= max_hits:
            break
    return hits


def _identity(feature: Feature) -> dict[str, str | None]:
    def first(key: str) -> str | None:
        values = feature.attributes.get(key)
        return values[0] if values else None

    gene_id = None
    for value in feature.attributes.get("Dbxref", ()):
        if value.startswith("GeneID:"):
            gene_id = value.split(":", 1)[1]
            break
    return {
        "feature_id": feature.feature_id,
        "name": first("Name"),
        "locus_tag": first("locus_tag"),
        "gene_id": gene_id,
    }


def map_feature_boundaries(
    feature: Feature,
    source_sequence: str,
    targets: dict[str, str],
    *,
    flank: int = 50,
) -> dict[str, Any]:
    start_center = feature.start
    end_center = feature.end - 1
    if start_center - flank < 0 or end_center + flank >= len(source_sequence):
        return {**_identity(feature), "mapping_status": "uncertain", "mapping_reason": "boundary anchor reaches source sequence edge"}
    start_anchor = source_sequence[start_center - flank : start_center + flank + 1]
    end_anchor = source_sequence[end_center - flank : end_center + flank + 1]
    start_hits = _hits(start_anchor, targets)
    end_hits = _hits(end_anchor, targets)
    base = {
        **_identity(feature),
        "source_seqid": feature.seqid,
        "source_start": feature.start,
        "source_end": feature.end,
        "source_strand": feature.strand,
        "anchor_length": flank * 2 + 1,
        "start_anchor_hit_count": len(start_hits),
        "end_anchor_hit_count": len(end_hits),
    }
    if not start_hits and not end_hits:
        return {**base, "mapping_status": "unmappable", "mapping_reason": "neither boundary anchor maps exactly"}
    if len(start_hits) != 1 or len(end_hits) != 1:
        return {**base, "mapping_status": "uncertain", "mapping_reason": "one or both boundary anchors are missing or non-unique"}
    left = start_hits[0]
    right = end_hits[0]
    if left["seqid"] != right["seqid"] or left["orientation"] != right["orientation"]:
        return {**base, "mapping_status": "conflict", "mapping_reason": "boundary anchors map to different target identities"}
    left_center = left["anchor_start"] + flank
    right_center = right["anchor_start"] + flank
    orientation = left["orientation"]
    if orientation == "+" and left_center <= right_center:
        mapped_start, mapped_end = left_center, right_center + 1
    elif orientation == "-" and right_center <= left_center:
        mapped_start, mapped_end = right_center, left_center + 1
    else:
        return {**base, "mapping_status": "conflict", "mapping_reason": "boundary anchor order conflicts with mapping orientation"}
    mapped_strand = feature.strand
    if orientation == "-" and feature.strand in {"+", "-"}:
        mapped_strand = "+" if feature.strand == "-" else "-"
    return {
        **base,
        "mapping_status": "consistent",
        "mapping_reason": "both 101-bp boundary anchors map uniquely and coherently",
        "target_seqid": left["seqid"],
        "target_start": mapped_start,
        "target_end": mapped_end,
        "assembly_orientation": orientation,
        "mapped_strand": mapped_strand,
        "start_boundary_offset": mapped_start - feature.start,
        "end_boundary_offset": mapped_end - feature.end,
    }


def select_probe_windows(features: list[Feature]) -> list[tuple[str, str, Feature]]:
    by_sequence: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        if feature.feature_type == "gene":
            by_sequence[feature.seqid].append(feature)
    selected: list[tuple[str, str, Feature]] = []
    for seqid, genes in sorted(by_sequence.items()):
        genes.sort(key=lambda item: (item.start, item.end, item.feature_id or ""))
        centers = {"first": 1, "middle": len(genes) // 2, "last": max(1, len(genes) - 2)}
        used: set[int] = set()
        for label, center in centers.items():
            indices = [max(0, center - 1), center, min(len(genes) - 1, center + 1)]
            window_id = f"{seqid}:{label}"
            for index in indices:
                if index not in used:
                    selected.append((window_id, label, genes[index]))
                    used.add(index)
    return selected


def run_mapping_probe(
    source_fasta_path: Path,
    source_gff_path: Path,
    target_fasta_path: Path,
    target_sequence_ids: set[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    source_sequences = read_fasta_sequences(source_fasta_path)
    target_sequences = {
        key: value
        for key, value in read_fasta_sequences(target_fasta_path).items()
        if key in target_sequence_ids
    }
    features, _, _ = read_gff3(source_gff_path)
    records: list[dict[str, Any]] = []
    for window_id, label, feature in select_probe_windows(features):
        record = map_feature_boundaries(feature, source_sequences[feature.seqid], target_sequences)
        record["window_id"] = window_id
        record["window_position"] = label
        records.append(record)

    by_window: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        by_window[record["window_id"]].append(record)
    for window_records in by_window.values():
        ordered = sorted(window_records, key=lambda item: item.get("source_start", -1))
        if not all(item["mapping_status"] == "consistent" for item in ordered):
            order_status = "uncertain"
        elif len({item["target_seqid"] for item in ordered}) != 1 or len({item["assembly_orientation"] for item in ordered}) != 1:
            order_status = "conflict"
        else:
            centers = [(item["target_start"] + item["target_end"]) // 2 for item in ordered]
            orientation = ordered[0]["assembly_orientation"]
            order_status = "preserved" if (centers == sorted(centers) if orientation == "+" else centers == sorted(centers, reverse=True)) else "conflict"
        for item in window_records:
            item["order_status"] = order_status
            item["adjacency_status"] = "unavailable"
            item["adjacency_reason"] = "exact boundary anchors establish collinear order only; target gene adjacency was not qualified"

    counts = Counter(record["mapping_status"] for record in records)
    for status in ("consistent", "conflict", "unmappable", "uncertain"):
        counts.setdefault(status, 0)
    source_sequences_covered = sorted({record.get("source_seqid") for record in records if record.get("source_seqid")})
    summary = {
        "probe_record_count": len(records),
        "mapping_status_counts": dict(sorted(counts.items())),
        "source_nuclear_sequences_covered": source_sequences_covered,
        "source_nuclear_sequence_count": len(source_sequences_covered),
        "window_count": len(by_window),
        "order_status_counts": dict(sorted(Counter(next(iter(items))["order_status"] for items in by_window.values()).items())),
        "adjacency_status_counts": dict(sorted(Counter(next(iter(items))["adjacency_status"] for items in by_window.values()).items())),
        "method": "exact 101-bp anchors centered on submitted gene start and end boundaries",
    }
    return sorted(records, key=lambda item: (item.get("source_seqid", ""), item.get("source_start", -1), item.get("feature_id") or "")), summary
