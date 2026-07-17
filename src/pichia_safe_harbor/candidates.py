from __future__ import annotations

from collections import defaultdict
from typing import Any, Callable, TypeVar

from .models import CandidateWindow, ExcludedRegion, ExclusionZone, IntergenicRegion, RuleFlag

RULE_VERSION = "slice1-placeholder-rules-v1"

T = TypeVar("T")


def _percentile(values: list[int], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    index = round((len(ordered) - 1) * fraction)
    return ordered[index]


def _count_by(items: list[T], key: Callable[[T], str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in items:
        counts[key(item)] = counts.get(key(item), 0) + 1
    return dict(sorted(counts.items()))


def _subtract_intervals(start: int, end: int, cuts: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Subtract zero or more [cut_start, cut_end) intervals from [start, end).

    Cuts may be unsorted or overlapping. Returns the remaining contiguous
    sub-intervals in ascending order; an empty list means the whole span was cut away.
    """
    if end <= start:
        return []
    if not cuts:
        return [(start, end)]
    merged: list[list[int]] = []
    for cut_start, cut_end in sorted(cuts):
        clipped_start = max(cut_start, start)
        clipped_end = min(cut_end, end)
        if clipped_end <= clipped_start:
            continue
        if merged and clipped_start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], clipped_end)
        else:
            merged.append([clipped_start, clipped_end])
    remaining: list[tuple[int, int]] = []
    cursor = start
    for cut_start, cut_end in merged:
        if cut_start > cursor:
            remaining.append((cursor, cut_start))
        cursor = max(cursor, cut_end)
    if cursor < end:
        remaining.append((cursor, end))
    return remaining


def _flags_for(region: IntergenicRegion, length: int, long_interval_threshold: float, long_interval_percentile: float) -> list[RuleFlag]:
    flags: list[RuleFlag] = []
    if region.orientation != "convergent":
        flags.append(
            RuleFlag(
                rule_id=f"neighbor_orientation_{region.orientation}",
                reason=f"neighboring gene orientation is {region.orientation}, not convergent",
            )
        )
    if long_interval_threshold > 0 and length > long_interval_threshold:
        flags.append(
            RuleFlag(
                rule_id="unusually_long_interval",
                reason=(
                    f"parent intergenic region length {length}bp is above the "
                    f"p{long_interval_percentile * 100:.0f} length threshold "
                    f"({long_interval_threshold}bp) and may contain an unannotated feature"
                ),
            )
        )
    return flags


def classify_candidate_windows(
    intergenic_regions: list[IntergenicRegion],
    buffer_distance_bp: int,
    min_candidate_window_bp: int,
    long_interval_percentile: float = 0.95,
    extra_exclusion_zones: list[ExclusionZone] = (),
) -> tuple[list[CandidateWindow], list[ExcludedRegion], dict[str, Any]]:
    """Apply the Slice 1 exclude/flag rules to a Slice 0 intergenic-region inventory.

    buffer_distance_bp and min_candidate_window_bp are caller-supplied parameters,
    not defaults: REQUIREMENTS.md 4.6 reserves threshold freezing for a dedicated ADR,
    so this function refuses to silently assume a value.

    extra_exclusion_zones (e.g. future repeat/centromere/telomere risk tracks) are
    optional; when one falls inside a region's buffered span, that region splits into
    as many independent candidate_window records as remain after subtraction, per
    REQUIREMENTS.md 4.3. With no zones (the default), behavior is unchanged from a
    single candidate window per surviving region.
    """
    if buffer_distance_bp < 0:
        raise ValueError("buffer_distance_bp must be >= 0")
    if min_candidate_window_bp < 1:
        raise ValueError("min_candidate_window_bp must be >= 1")

    long_interval_threshold = _percentile(
        [region.length for region in intergenic_regions], long_interval_percentile
    )

    zones_by_seqid: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    for zone in extra_exclusion_zones:
        zones_by_seqid[zone.seqid].append((zone.start, zone.end, zone.rule_id, zone.reason))

    windows: list[CandidateWindow] = []
    excluded: list[ExcludedRegion] = []
    for region in sorted(intergenic_regions, key=lambda item: (item.seqid, item.start, item.region_id)):
        if region.left_boundary_confidence != "high" or region.right_boundary_confidence != "high":
            excluded.append(
                ExcludedRegion(
                    region_id=f"{region.region_id}:excluded",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=region.start,
                    end=region.end,
                    length=region.length,
                    rule_id="boundary_confidence_insufficient",
                    reason=(
                        f"left_boundary_confidence={region.left_boundary_confidence}, "
                        f"right_boundary_confidence={region.right_boundary_confidence}; the "
                        "conservative candidate set requires high confidence on both sides "
                        "pending independent annotation support (ADR-0004/ADR-0005)"
                    ),
                )
            )
            continue

        candidate_start = region.start + buffer_distance_bp
        candidate_end = region.end - buffer_distance_bp

        if candidate_end <= candidate_start:
            excluded.append(
                ExcludedRegion(
                    region_id=f"{region.region_id}:excluded",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=region.start,
                    end=region.end,
                    length=region.length,
                    rule_id="below_minimum_window_length",
                    reason=(
                        f"remaining {candidate_end - candidate_start}bp after {buffer_distance_bp}bp "
                        f"buffer on each side < minimum {min_candidate_window_bp}bp"
                    ),
                )
            )
            continue

        cuts = [
            (start, end, rule_id, reason)
            for start, end, rule_id, reason in zones_by_seqid.get(region.seqid, [])
            if end > candidate_start and start < candidate_end
        ]

        if not cuts:
            remaining = candidate_end - candidate_start
            if remaining < min_candidate_window_bp:
                excluded.append(
                    ExcludedRegion(
                        region_id=f"{region.region_id}:excluded",
                        parent_region_id=region.region_id,
                        seqid=region.seqid,
                        start=region.start,
                        end=region.end,
                        length=region.length,
                        rule_id="below_minimum_window_length",
                        reason=(
                            f"remaining {remaining}bp after {buffer_distance_bp}bp buffer on each side "
                            f"< minimum {min_candidate_window_bp}bp"
                        ),
                    )
                )
                continue
            flags = _flags_for(region, remaining, long_interval_threshold, long_interval_percentile)
            windows.append(
                CandidateWindow(
                    candidate_id=f"{region.region_id}:window:001",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=candidate_start,
                    end=candidate_end,
                    length=remaining,
                    left_entity_ids=region.left_entity_ids,
                    right_entity_ids=region.right_entity_ids,
                    orientation=region.orientation,
                    left_boundary_confidence=region.left_boundary_confidence,
                    right_boundary_confidence=region.right_boundary_confidence,
                    evidence_level="predicted",
                    structural_tier="convergent_clean" if not flags else "flagged",
                    split_index=1,
                    split_count=1,
                    rule_flags=tuple(flags),
                )
            )
            continue

        for start, end, rule_id, reason in cuts:
            clipped_start, clipped_end = max(start, candidate_start), min(end, candidate_end)
            excluded.append(
                ExcludedRegion(
                    region_id=f"{region.region_id}:excluded:{rule_id}:{clipped_start}",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=clipped_start,
                    end=clipped_end,
                    length=clipped_end - clipped_start,
                    rule_id=rule_id,
                    reason=reason,
                )
            )

        remaining_spans = _subtract_intervals(
            candidate_start, candidate_end, [(start, end) for start, end, _, _ in cuts]
        )
        valid_spans = [span for span in remaining_spans if span[1] - span[0] >= min_candidate_window_bp]
        short_spans = [span for span in remaining_spans if span[1] - span[0] < min_candidate_window_bp]
        for span_start, span_end in short_spans:
            excluded.append(
                ExcludedRegion(
                    region_id=f"{region.region_id}:excluded:below_minimum_window_length:{span_start}",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=span_start,
                    end=span_end,
                    length=span_end - span_start,
                    rule_id="below_minimum_window_length",
                    reason=(
                        f"remaining fragment {span_end - span_start}bp after buffer and exclusion-zone "
                        f"subtraction < minimum {min_candidate_window_bp}bp"
                    ),
                )
            )

        split_count = len(valid_spans)
        for index, (span_start, span_end) in enumerate(valid_spans, start=1):
            span_length = span_end - span_start
            flags = _flags_for(region, span_length, long_interval_threshold, long_interval_percentile)
            windows.append(
                CandidateWindow(
                    candidate_id=f"{region.region_id}:window:{index:03d}",
                    parent_region_id=region.region_id,
                    seqid=region.seqid,
                    start=span_start,
                    end=span_end,
                    length=span_length,
                    left_entity_ids=region.left_entity_ids if span_start == candidate_start else (),
                    right_entity_ids=region.right_entity_ids if span_end == candidate_end else (),
                    orientation=region.orientation,
                    left_boundary_confidence=region.left_boundary_confidence,
                    right_boundary_confidence=region.right_boundary_confidence,
                    evidence_level="predicted",
                    structural_tier="convergent_clean" if not flags else "flagged",
                    split_index=index,
                    split_count=split_count,
                    rule_flags=tuple(flags),
                )
            )

    summary = {
        "input_intergenic_region_count": len(intergenic_regions),
        "candidate_window_count": len(windows),
        "excluded_region_count": len(excluded),
        "excluded_by_rule": _count_by(excluded, lambda item: item.rule_id),
        "candidate_structural_tier_counts": _count_by(windows, lambda item: item.structural_tier),
        "long_interval_length_threshold_bp": long_interval_threshold,
    }
    return windows, excluded, summary
