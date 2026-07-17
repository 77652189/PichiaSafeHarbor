from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from .errors import ContractError
from .models import (
    BoundaryCluster,
    Diagnostics,
    Feature,
    FunctionalEntity,
    IntergenicRegion,
    SequenceInfo,
    TerminalRegion,
)

BOUNDARY_TYPES = {
    "gene",
    "pseudogene",
    "ncRNA_gene",
    "tRNA_gene",
    "rRNA_gene",
    "snRNA_gene",
    "snoRNA_gene",
}
TRANSCRIPT_TYPES = {
    "mRNA",
    "transcript",
    "ncRNA",
    "lnc_RNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
}
CHILD_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}


def map_annotation_sequences(
    features: list[Feature],
    fasta: dict[str, SequenceInfo],
    explicit_map: dict[str, str],
) -> tuple[list[Feature], dict[str, str]]:
    annotation_seqids = sorted({feature.seqid for feature in features})
    mapping: dict[str, str] = {}
    reverse: dict[str, str] = {}
    for seqid in annotation_seqids:
        target = explicit_map.get(seqid, seqid if seqid in fasta else "")
        if not target or target not in fasta:
            raise ContractError(f"annotation sequence {seqid!r} cannot map to FASTA uniquely")
        if target in reverse and reverse[target] != seqid:
            raise ContractError(
                f"annotation sequences {reverse[target]!r} and {seqid!r} both map to FASTA {target!r}"
            )
        mapping[seqid] = target
        reverse[target] = seqid
    mapped = [
        Feature(
            seqid=mapping[feature.seqid],
            source=feature.source,
            feature_type=feature.feature_type,
            start=feature.start,
            end=feature.end,
            strand=feature.strand,
            attributes=feature.attributes,
            line_number=feature.line_number,
        )
        for feature in features
    ]
    return mapped, mapping


def validate_coordinates(features: list[Feature], fasta: dict[str, SequenceInfo]) -> None:
    for feature in features:
        if feature.seqid not in fasta:
            raise ContractError(f"feature sequence is absent from FASTA: {feature.seqid}")
        if feature.start < 0 or feature.end > fasta[feature.seqid].length:
            raise ContractError(
                f"annotation coordinate out of bounds at line {feature.line_number}: "
                f"{feature.seqid}:{feature.start}-{feature.end} / {fasta[feature.seqid].length}"
            )


def validate_identity(reference: dict[str, Any], directives: dict[str, str]) -> None:
    expected = reference["assembly_accession"]
    declared = " ".join(directives.values())
    if reference.get("require_gff_assembly_directive", True) and expected not in declared:
        raise ContractError(
            f"annotation does not declare expected assembly {expected}; directives={directives}"
        )


def normalize_entities(features: list[Feature]) -> tuple[list[FunctionalEntity], Diagnostics]:
    diagnostics = Diagnostics()
    id_records: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        if feature.feature_id:
            id_records[feature.feature_id].append(feature)
    for feature_id, records in sorted(id_records.items()):
        if len(records) > 1:
            diagnostics.duplicate_ids.append(
                {"id": feature_id, "lines": [record.line_number for record in records]}
            )
            signatures = {
                (
                    record.seqid,
                    record.feature_type,
                    record.start,
                    record.end,
                    record.strand,
                    record.parents,
                )
                for record in records
            }
            compatible_segmented_child = (
                all(record.feature_type in CHILD_TYPES for record in records)
                and len(
                    {
                        (record.seqid, record.feature_type, record.strand, record.parents)
                        for record in records
                    }
                )
                == 1
            )
            if len(signatures) > 1 and not compatible_segmented_child:
                raise ContractError(f"conflicting duplicate annotation ID: {feature_id}")

    known_ids = set(id_records)
    for feature in features:
        if feature.feature_type in BOUNDARY_TYPES | TRANSCRIPT_TYPES and not feature.feature_id:
            diagnostics.missing_ids.append(
                {"line": feature.line_number, "type": feature.feature_type, "seqid": feature.seqid}
            )
        for parent in feature.parents:
            if parent not in known_ids:
                diagnostics.missing_parents.append(
                    {"line": feature.line_number, "feature_id": feature.feature_id, "parent": parent}
                )
            elif len(id_records[parent]) > 1:
                diagnostics.ambiguous_parents.append(
                    {"line": feature.line_number, "feature_id": feature.feature_id, "parent": parent}
                )

    entities: list[FunctionalEntity] = []
    boundary_ids = {
        feature.feature_id
        for feature in features
        if feature.feature_type in BOUNDARY_TYPES and feature.feature_id
    }
    parent_graph = {
        feature.feature_id: feature.parents
        for feature in features
        if feature.feature_id and len(id_records[feature.feature_id]) == 1
    }

    def has_ancestor(feature: Feature, ancestors: set[str]) -> bool:
        pending = list(feature.parents)
        visited: set[str] = set()
        while pending:
            parent = pending.pop()
            if parent in ancestors:
                return True
            if parent in visited:
                continue
            visited.add(parent)
            pending.extend(parent_graph.get(parent, ()))
        return False

    standalone_transcript_ids = {
        feature.feature_id
        for feature in features
        if feature.feature_type in TRANSCRIPT_TYPES
        and feature.feature_id
        and not has_ancestor(feature, boundary_ids)
    }
    emitted_exact_ids: set[str] = set()

    def attribute_values(records: list[Feature], key: str) -> tuple[str, ...]:
        return tuple(
            dict.fromkeys(
                value
                for record in records
                for value in record.attributes.get(key, ())
            )
        )

    def is_partial(records: list[Feature]) -> bool:
        return any(
            value.lower() == "true"
            for record in records
            for value in record.attributes.get("partial", ())
        )

    def confidence(records: list[Feature], range_key: str) -> str:
        if attribute_values(records, range_key):
            return "uncertain"
        has_any_range = bool(
            attribute_values(records, "start_range")
            or attribute_values(records, "end_range")
        )
        if is_partial(records) and not has_any_range:
            return "uncertain"
        return "high"

    for index, feature in enumerate(features, 1):
        if feature.feature_type in BOUNDARY_TYPES:
            base_id = feature.feature_id or f"missing-id-line-{feature.line_number}"
        elif feature.feature_type in TRANSCRIPT_TYPES:
            if has_ancestor(feature, boundary_ids):
                continue
            base_id = feature.feature_id or f"orphan-transcript-line-{feature.line_number}"
        elif feature.feature_type in CHILD_TYPES:
            if has_ancestor(feature, boundary_ids | standalone_transcript_ids):
                continue
            base_id = feature.feature_id or f"orphan-{feature.feature_type}-line-{feature.line_number}"
        else:
            diagnostics.ignored_feature_types[feature.feature_type] = (
                diagnostics.ignored_feature_types.get(feature.feature_type, 0) + 1
            )
            continue
        if feature.feature_id and feature.feature_id in emitted_exact_ids:
            continue
        if feature.feature_id:
            emitted_exact_ids.add(feature.feature_id)
        duplicate_suffix = ""
        span_records = id_records.get(feature.feature_id, [feature]) if feature.feature_id else [feature]
        span_start = min(record.start for record in span_records)
        span_end = max(record.end for record in span_records)
        partial = is_partial(span_records)
        start_range = attribute_values(span_records, "start_range")
        end_range = attribute_values(span_records, "end_range")
        start_confidence = confidence(span_records, "start_range")
        end_confidence = confidence(span_records, "end_range")
        if feature.strand == "-":
            five_prime_confidence = end_confidence
            three_prime_confidence = start_confidence
        else:
            five_prime_confidence = start_confidence
            three_prime_confidence = end_confidence
        entities.append(
            FunctionalEntity(
                entity_id=base_id + duplicate_suffix,
                seqid=feature.seqid,
                start=span_start,
                end=span_end,
                strand=feature.strand,
                entity_type=feature.feature_type,
                source_feature_ids=(feature.feature_id,) if feature.feature_id else (),
                partial=partial,
                start_range=start_range,
                end_range=end_range,
                genomic_start_confidence=start_confidence,
                genomic_end_confidence=end_confidence,
                five_prime_confidence=five_prime_confidence,
                three_prime_confidence=three_prime_confidence,
            )
        )

    by_seqid: dict[str, list[FunctionalEntity]] = defaultdict(list)
    for entity in entities:
        by_seqid[entity.seqid].append(entity)
    for seqid, sequence_entities in by_seqid.items():
        original_order = [(entity.start, entity.end, entity.entity_id) for entity in sequence_entities]
        sorted_order = sorted(original_order)
        if original_order != sorted_order:
            diagnostics.abnormal_order.append({"seqid": seqid, "entity_count": len(sequence_entities)})
        ordered = sorted(sequence_entities, key=lambda item: (item.start, item.end, item.entity_id))
        for left_index, left in enumerate(ordered):
            for right in ordered[left_index + 1 :]:
                if right.start >= left.end:
                    break
                kind = "nested" if right.end <= left.end else "overlaps"
                getattr(diagnostics, kind).append(
                    {
                        "seqid": seqid,
                        "left": left.entity_id,
                        "right": right.entity_id,
                        "intersection": [right.start, min(left.end, right.end)],
                    }
                )
    return sorted(entities, key=lambda item: (item.seqid, item.start, item.end, item.entity_id)), diagnostics


def merge_boundaries(entities: list[FunctionalEntity]) -> list[BoundaryCluster]:
    if not entities:
        return []
    ordered = sorted(entities, key=lambda item: (item.start, item.end, item.entity_id))
    clusters: list[BoundaryCluster] = []
    start = ordered[0].start
    end = ordered[0].end
    members = [ordered[0]]
    seqid = ordered[0].seqid

    def make_cluster(cluster_members: list[FunctionalEntity]) -> BoundaryCluster:
        cluster_start = min(item.start for item in cluster_members)
        cluster_end = max(item.end for item in cluster_members)
        left_members = [item for item in cluster_members if item.start == cluster_start]
        right_members = [item for item in cluster_members if item.end == cluster_end]
        return BoundaryCluster(
            seqid=seqid,
            start=cluster_start,
            end=cluster_end,
            entity_ids=tuple(item.entity_id for item in cluster_members),
            left_boundary_entity_ids=tuple(item.entity_id for item in left_members),
            right_boundary_entity_ids=tuple(item.entity_id for item in right_members),
            left_boundary_strands=tuple(item.strand for item in left_members),
            right_boundary_strands=tuple(item.strand for item in right_members),
            left_boundary_confidences=tuple(
                item.genomic_start_confidence for item in left_members
            ),
            right_boundary_confidences=tuple(
                item.genomic_end_confidence for item in right_members
            ),
        )

    for entity in ordered[1:]:
        if entity.seqid != seqid:
            raise ContractError("merge_boundaries received entities from multiple sequences")
        if entity.start < end:
            end = max(end, entity.end)
            members.append(entity)
        else:
            clusters.append(make_cluster(members))
            start, end = entity.start, entity.end
            members = [entity]
    clusters.append(make_cluster(members))
    return clusters


def orientation(left: str, right: str) -> str:
    if left not in {"+", "-"} or right not in {"+", "-"}:
        return "unknown"
    if left == right:
        return "tandem"
    if left == "+" and right == "-":
        return "convergent"
    return "divergent"


def build_inventory(
    sequences: dict[str, SequenceInfo],
    entities: list[FunctionalEntity],
    sequence_classes: dict[str, str],
    sequence_endpoints: dict[str, dict[str, Any]],
) -> tuple[list[IntergenicRegion], list[TerminalRegion], dict[str, Any]]:
    unknown = sorted(set(sequences) - set(sequence_classes))
    extra = sorted(set(sequence_classes) - set(sequences))
    if unknown or extra:
        raise ContractError(f"sequence classification mismatch: unclassified={unknown}, absent={extra}")
    endpoint_unknown = sorted(set(sequences) - set(sequence_endpoints))
    endpoint_extra = sorted(set(sequence_endpoints) - set(sequences))
    if endpoint_unknown or endpoint_extra:
        raise ContractError(
            "sequence endpoint completeness mismatch: "
            f"unclassified={endpoint_unknown}, absent={endpoint_extra}"
        )
    by_sequence: dict[str, list[FunctionalEntity]] = defaultdict(list)
    for entity in entities:
        by_sequence[entity.seqid].append(entity)

    intergenic: list[IntergenicRegion] = []
    terminal: list[TerminalRegion] = []
    sequence_stats: dict[str, Any] = {}
    for seqid in sorted(sequences):
        seq_class = sequence_classes[seqid]
        endpoint_record = sequence_endpoints[seqid]
        five_prime = endpoint_record["five_prime"]
        three_prime = endpoint_record["three_prime"]
        for endpoint_name, endpoint in (("five_prime", five_prime), ("three_prime", three_prime)):
            if endpoint.get("status") not in {"complete", "partial", "unknown"}:
                raise ContractError(
                    f"invalid {endpoint_name} completeness for {seqid}: {endpoint.get('status')}"
                )
        sequence_entities = by_sequence.get(seqid, [])
        if seq_class != "nuclear_chromosome":
            sequence_stats[seqid] = {
                "sequence_class": seq_class,
                "length": sequences[seqid].length,
                "functional_entity_count": len(sequence_entities),
                "boundary_cluster_count": 0,
                "intergenic_region_count": 0,
                "excluded_from_nuclear_statistics": True,
                "endpoint_completeness": endpoint_record,
            }
            continue
        clusters = merge_boundaries(sequence_entities)
        if not clusters:
            statuses = {five_prime["status"], three_prime["status"]}
            combined_status = next(iter(statuses)) if len(statuses) == 1 else "unknown"
            terminal.append(
                TerminalRegion(
                    region_id=f"{seqid}:whole-sequence",
                    seqid=seqid,
                    start=0,
                    end=sequences[seqid].length,
                    length=sequences[seqid].length,
                    side="whole_sequence_no_features",
                    five_prime_completeness=five_prime["status"],
                    three_prime_completeness=three_prime["status"],
                    adjacent_sequence_end="both",
                    adjacent_end_completeness=combined_status,
                    adjacent_end_evidence_source=(
                        f"5p: {five_prime['evidence_source']} | 3p: {three_prime['evidence_source']}"
                    ),
                    adjacent_end_original_declaration=(
                        f"5p: {five_prime['original_declaration']} | "
                        f"3p: {three_prime['original_declaration']}"
                    ),
                    interpretable_as_true_telomere_distance=(
                        five_prime["status"] == "complete"
                        and three_prime["status"] == "complete"
                    ),
                )
            )
        else:
            first = clusters[0]
            if first.start > 0:
                terminal.append(
                    TerminalRegion(
                        region_id=f"{seqid}:5p-terminal",
                        seqid=seqid,
                        start=0,
                        end=first.start,
                        length=first.start,
                        side="5p",
                        neighbor_entity_ids=first.left_boundary_entity_ids,
                        five_prime_completeness=five_prime["status"],
                        three_prime_completeness=three_prime["status"],
                        adjacent_sequence_end="5p",
                        adjacent_end_completeness=five_prime["status"],
                        adjacent_end_evidence_source=five_prime["evidence_source"],
                        adjacent_end_original_declaration=five_prime["original_declaration"],
                        interpretable_as_true_telomere_distance=five_prime["status"] == "complete",
                    )
                )
            for index, (left, right) in enumerate(zip(clusters, clusters[1:]), 1):
                if right.start <= left.end:
                    continue
                gap = right.start - left.end
                intergenic.append(
                    IntergenicRegion(
                        region_id=f"{seqid}:intergenic:{index:06d}",
                        seqid=seqid,
                        start=left.end,
                        end=right.start,
                        length=gap,
                        left_entity_ids=left.right_boundary_entity_ids,
                        right_entity_ids=right.left_boundary_entity_ids,
                        left_strand=left.right_boundary_strand,
                        right_strand=right.left_boundary_strand,
                        orientation=orientation(left.right_boundary_strand, right.left_boundary_strand),
                        gap_length_bp=gap,
                        left_boundary_confidence=left.right_boundary_confidence,
                        right_boundary_confidence=right.left_boundary_confidence,
                    )
                )
            last = clusters[-1]
            if last.end < sequences[seqid].length:
                terminal.append(
                    TerminalRegion(
                        region_id=f"{seqid}:3p-terminal",
                        seqid=seqid,
                        start=last.end,
                        end=sequences[seqid].length,
                        length=sequences[seqid].length - last.end,
                        side="3p",
                        neighbor_entity_ids=last.right_boundary_entity_ids,
                        five_prime_completeness=five_prime["status"],
                        three_prime_completeness=three_prime["status"],
                        adjacent_sequence_end="3p",
                        adjacent_end_completeness=three_prime["status"],
                        adjacent_end_evidence_source=three_prime["evidence_source"],
                        adjacent_end_original_declaration=three_prime["original_declaration"],
                        interpretable_as_true_telomere_distance=three_prime["status"] == "complete",
                    )
                )
        sequence_stats[seqid] = {
            "sequence_class": seq_class,
            "length": sequences[seqid].length,
            "functional_entity_count": len(sequence_entities),
            "boundary_cluster_count": len(clusters),
            "intergenic_region_count": sum(1 for item in intergenic if item.seqid == seqid),
            "excluded_from_nuclear_statistics": False,
            "endpoint_completeness": endpoint_record,
        }

    orientation_counts = Counter(item.orientation for item in intergenic)
    for category in ("convergent", "tandem", "divergent", "unknown"):
        orientation_counts.setdefault(category, 0)
    orientation_lengths: dict[str, list[int]] = defaultdict(list)
    for item in intergenic:
        orientation_lengths[item.orientation].append(item.length)
    high_confidence = [
        item
        for item in intergenic
        if item.left_boundary_confidence == "high"
        and item.right_boundary_confidence == "high"
    ]
    high_orientation_counts = Counter(item.orientation for item in high_confidence)
    for category in ("convergent", "tandem", "divergent", "unknown"):
        high_orientation_counts.setdefault(category, 0)
    partial_entities = [item for item in entities if item.partial]
    start_range_entities = [item for item in entities if item.start_range]
    end_range_entities = [item for item in entities if item.end_range]
    entity_type_boundary_summary: dict[str, dict[str, int]] = {}
    for entity_type in sorted({item.entity_type for item in entities}):
        typed = [item for item in entities if item.entity_type == entity_type]
        entity_type_boundary_summary[entity_type] = {
            "count": len(typed),
            "partial_count": sum(1 for item in typed if item.partial),
            "start_range_count": sum(1 for item in typed if item.start_range),
            "end_range_count": sum(1 for item in typed if item.end_range),
        }
    summary = {
        "sequence_statistics": sequence_stats,
        "nuclear_chromosome_count": sum(
            1 for value in sequence_classes.values() if value == "nuclear_chromosome"
        ),
        "intergenic_region_count": len(intergenic),
        "orientation_counts": dict(sorted(orientation_counts.items())),
        "intergenic_length_distribution": describe([item.length for item in intergenic]),
        "orientation_length_distributions": {
            key: describe(orientation_lengths.get(key, []))
            for key in ("convergent", "tandem", "divergent", "unknown")
        },
        "high_confidence_intergenic": {
            "count": len(high_confidence),
            "fraction": round(len(high_confidence) / len(intergenic), 6) if intergenic else 0.0,
            "orientation_counts": dict(sorted(high_orientation_counts.items())),
            "length_distribution": describe([item.length for item in high_confidence]),
        },
        "annotation_boundary_summary": {
            "functional_entity_count": len(entities),
            "partial_entity_count": len(partial_entities),
            "start_range_entity_count": len(start_range_entities),
            "end_range_entity_count": len(end_range_entities),
            "partial_fraction": round(len(partial_entities) / len(entities), 6) if entities else 0.0,
            "by_entity_type": entity_type_boundary_summary,
        },
    }
    return intergenic, terminal, summary


def describe(values: list[int]) -> dict[str, int | float | None]:
    if not values:
        return {"count": 0, "min": None, "p25": None, "median": None, "p75": None, "max": None, "mean": None}
    ordered = sorted(values)

    def percentile(fraction: float) -> int:
        index = round((len(ordered) - 1) * fraction)
        return ordered[index]

    return {
        "count": len(ordered),
        "min": ordered[0],
        "p25": percentile(0.25),
        "median": percentile(0.5),
        "p75": percentile(0.75),
        "max": ordered[-1],
        "mean": round(sum(ordered) / len(ordered), 6),
    }
