from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class SequenceInfo:
    name: str
    length: int
    description: str = ""


@dataclass(frozen=True)
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict[str, tuple[str, ...]]
    line_number: int

    @property
    def feature_id(self) -> str | None:
        values = self.attributes.get("ID")
        return values[0] if values else None

    @property
    def parents(self) -> tuple[str, ...]:
        return self.attributes.get("Parent", ())


@dataclass(frozen=True)
class FunctionalEntity:
    entity_id: str
    seqid: str
    start: int
    end: int
    strand: str
    entity_type: str
    source_feature_ids: tuple[str, ...] = ()
    partial: bool = False
    start_range: tuple[str, ...] = ()
    end_range: tuple[str, ...] = ()
    genomic_start_confidence: str = "high"
    genomic_end_confidence: str = "high"
    five_prime_confidence: str = "high"
    three_prime_confidence: str = "high"

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["source_feature_ids"] = list(self.source_feature_ids)
        value["start_range"] = list(self.start_range)
        value["end_range"] = list(self.end_range)
        return value


@dataclass(frozen=True)
class BoundaryCluster:
    seqid: str
    start: int
    end: int
    entity_ids: tuple[str, ...]
    left_boundary_entity_ids: tuple[str, ...]
    right_boundary_entity_ids: tuple[str, ...]
    left_boundary_strands: tuple[str, ...]
    right_boundary_strands: tuple[str, ...]
    left_boundary_confidences: tuple[str, ...]
    right_boundary_confidences: tuple[str, ...]

    @staticmethod
    def _representative_strand(strands: tuple[str, ...]) -> str:
        known = {strand for strand in strands if strand in {"+", "-"}}
        return next(iter(known)) if len(known) == 1 else "."

    @property
    def left_boundary_strand(self) -> str:
        return self._representative_strand(self.left_boundary_strands)

    @property
    def right_boundary_strand(self) -> str:
        return self._representative_strand(self.right_boundary_strands)

    @staticmethod
    def _representative_confidence(confidences: tuple[str, ...]) -> str:
        if confidences and all(value == "high" for value in confidences):
            return "high"
        if any(value == "uncertain" for value in confidences):
            return "uncertain"
        return "unknown"

    @property
    def left_boundary_confidence(self) -> str:
        return self._representative_confidence(self.left_boundary_confidences)

    @property
    def right_boundary_confidence(self) -> str:
        return self._representative_confidence(self.right_boundary_confidences)


@dataclass(frozen=True)
class IntergenicRegion:
    region_id: str
    seqid: str
    start: int
    end: int
    length: int
    left_entity_ids: tuple[str, ...]
    right_entity_ids: tuple[str, ...]
    left_strand: str
    right_strand: str
    orientation: str
    gap_length_bp: int
    left_boundary_confidence: str
    right_boundary_confidence: str
    coordinate_system: str = "0-based-half-open"

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["left_entity_ids"] = list(self.left_entity_ids)
        value["right_entity_ids"] = list(self.right_entity_ids)
        return value


@dataclass(frozen=True)
class TerminalRegion:
    region_id: str
    seqid: str
    start: int
    end: int
    length: int
    side: str
    neighbor_entity_ids: tuple[str, ...] = ()
    five_prime_completeness: str = "unknown"
    three_prime_completeness: str = "unknown"
    adjacent_sequence_end: str = "unknown"
    adjacent_end_completeness: str = "unknown"
    adjacent_end_evidence_source: str = "unavailable"
    adjacent_end_original_declaration: str = "unavailable"
    interpretable_as_true_telomere_distance: bool = False
    coordinate_system: str = "0-based-half-open"

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["neighbor_entity_ids"] = list(self.neighbor_entity_ids)
        return value


@dataclass(frozen=True)
class RuleFlag:
    rule_id: str
    reason: str


@dataclass(frozen=True)
class ExclusionZone:
    seqid: str
    start: int
    end: int
    rule_id: str
    reason: str


@dataclass(frozen=True)
class CandidateWindow:
    candidate_id: str
    parent_region_id: str
    seqid: str
    start: int
    end: int
    length: int
    left_entity_ids: tuple[str, ...]
    right_entity_ids: tuple[str, ...]
    orientation: str
    left_boundary_confidence: str
    right_boundary_confidence: str
    evidence_level: str
    structural_tier: str
    sequence: str = ""
    split_index: int = 1
    split_count: int = 1
    rule_flags: tuple[RuleFlag, ...] = ()
    coordinate_system: str = "0-based-half-open"

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["left_entity_ids"] = list(self.left_entity_ids)
        value["right_entity_ids"] = list(self.right_entity_ids)
        value["rule_flags"] = [asdict(flag) for flag in self.rule_flags]
        return value


@dataclass(frozen=True)
class ExcludedRegion:
    region_id: str
    parent_region_id: str
    seqid: str
    start: int
    end: int
    length: int
    rule_id: str
    reason: str
    coordinate_system: str = "0-based-half-open"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Diagnostics:
    duplicate_ids: list[dict[str, Any]] = field(default_factory=list)
    missing_ids: list[dict[str, Any]] = field(default_factory=list)
    missing_parents: list[dict[str, Any]] = field(default_factory=list)
    ambiguous_parents: list[dict[str, Any]] = field(default_factory=list)
    abnormal_order: list[dict[str, Any]] = field(default_factory=list)
    overlaps: list[dict[str, Any]] = field(default_factory=list)
    nested: list[dict[str, Any]] = field(default_factory=list)
    ignored_feature_types: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
