from __future__ import annotations

from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.parsers import read_fasta_sequences, read_gff3


def test_gff_closed_to_half_open_conversion(tmp_path: Path) -> None:
    path = tmp_path / "coords.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr\tx\tgene\t1\t1\t.\t+\t.\tID=first\n"
        "chr\tx\tgene\t10\t10\t.\t-\t.\tID=last\n"
        "chr\tx\tgene\t1\t10\t.\t+\t.\tID=all\n",
        encoding="utf-8",
    )
    features, _, _ = read_gff3(path)
    assert [(item.start, item.end) for item in features] == [(0, 1), (9, 10), (0, 10)]
    assert [item.end - item.start for item in features] == [1, 1, 10]


@pytest.mark.parametrize("coords", ["0\t1", "-1\t1", "3\t2"])
def test_invalid_gff_coordinates_fail(tmp_path: Path, coords: str) -> None:
    path = tmp_path / "bad.gff3"
    path.write_text(f"chr\tx\tgene\t{coords}\t.\t+\t.\tID=bad\n", encoding="utf-8")
    with pytest.raises(ContractError):
        read_gff3(path)


def test_read_fasta_sequences_joins_wrapped_lines(tmp_path: Path) -> None:
    path = tmp_path / "seqs.fna"
    path.write_text(">chr1 first sequence\nACGT\nACGT\n>chr2 second\nTTTT\n", encoding="utf-8")
    sequences = read_fasta_sequences(path)
    assert sequences == {"chr1": "ACGTACGT", "chr2": "TTTT"}


def test_read_fasta_sequences_rejects_duplicate_name(tmp_path: Path) -> None:
    path = tmp_path / "dup.fna"
    path.write_text(">chr1\nACGT\n>chr1\nTTTT\n", encoding="utf-8")
    with pytest.raises(ContractError, match="duplicate FASTA sequence name"):
        read_fasta_sequences(path)


def test_read_fasta_sequences_rejects_invalid_symbol(tmp_path: Path) -> None:
    path = tmp_path / "bad.fna"
    path.write_text(">chr1\nACGTXQ\n", encoding="utf-8")
    with pytest.raises(ContractError, match="invalid FASTA symbols"):
        read_fasta_sequences(path)

