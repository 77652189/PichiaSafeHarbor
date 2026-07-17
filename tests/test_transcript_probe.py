from __future__ import annotations

from pichia_safe_harbor.transcript_probe import _parse_sam, _reference_identity


def test_parse_sam_preserves_mapping_and_multiplicity() -> None:
    record = _parse_sam(
        "read1\t419\tchr1\t5\t0\t4M2N4M\t=\t20\t15\tACGTACGT\t????????\tNH:i:2\tNM:i:1"
    )
    assert record["read_id"] == "read1"
    assert record["source_position"] == 5
    assert record["nh"] == 2
    assert record["nm"] == 1


def test_reference_identity_handles_splice_cigar() -> None:
    record = {
        "source_position": 1,
        "cigar": "4M2N4M",
        "read_sequence": "ACGTGTAC",
    }
    result = _reference_identity(record, "ACGTNNGTAC")
    assert result["coordinate_valid"] is True
    assert result["reference_identity"] == 1.0
    assert result["splice_junctions"] == 1
    assert result["target_end"] == 10


def test_reference_identity_rejects_out_of_bounds_alignment() -> None:
    record = {"source_position": 5, "cigar": "8M", "read_sequence": "ACGTACGT"}
    result = _reference_identity(record, "AAAAAA")
    assert result["coordinate_valid"] is False
    assert result["reference_identity"] is None
