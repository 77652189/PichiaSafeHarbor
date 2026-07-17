from __future__ import annotations

from collections import Counter
from pathlib import Path
from urllib.parse import unquote

from .errors import ContractError
from .io_utils import open_text
from .models import Feature, SequenceInfo


def read_fasta_index(path: Path) -> dict[str, SequenceInfo]:
    sequences: dict[str, SequenceInfo] = {}
    current_name: str | None = None
    current_description = ""
    current_length = 0
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = SequenceInfo(
                        current_name, current_length, current_description
                    )
                header = line[1:].strip()
                current_name = header.split(None, 1)[0]
                if current_name in sequences:
                    raise ContractError(f"duplicate FASTA sequence name: {current_name}")
                current_description = header
                current_length = 0
                continue
            if current_name is None:
                raise ContractError(f"FASTA sequence appears before header at line {line_number}")
            sequence = "".join(line.split()).upper()
            invalid = set(sequence) - set("ACGTURYSWKMBDHVN.-*")
            if invalid:
                raise ContractError(
                    f"invalid FASTA symbols at line {line_number}: {''.join(sorted(invalid))}"
                )
            current_length += len(sequence)
    if current_name is not None:
        sequences[current_name] = SequenceInfo(current_name, current_length, current_description)
    if not sequences:
        raise ContractError("FASTA contains no sequences")
    return sequences


def read_fasta_sequences(path: Path) -> dict[str, str]:
    """Read full sequence content per seqid (unlike read_fasta_index, which only
    computes lengths). Kept independent of read_fasta_index rather than sharing
    its loop, to avoid any risk of changing the already-tested baseline path."""
    sequences: dict[str, list[str]] = {}
    current_name: str | None = None
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                current_name = header.split(None, 1)[0]
                if current_name in sequences:
                    raise ContractError(f"duplicate FASTA sequence name: {current_name}")
                sequences[current_name] = []
                continue
            if current_name is None:
                raise ContractError(f"FASTA sequence appears before header at line {line_number}")
            segment = "".join(line.split()).upper()
            invalid = set(segment) - set("ACGTURYSWKMBDHVN.-*")
            if invalid:
                raise ContractError(
                    f"invalid FASTA symbols at line {line_number}: {''.join(sorted(invalid))}"
                )
            sequences[current_name].append(segment)
    if not sequences:
        raise ContractError("FASTA contains no sequences")
    return {name: "".join(parts) for name, parts in sequences.items()}


def _parse_attributes(raw: str) -> dict[str, tuple[str, ...]]:
    attributes: dict[str, list[str]] = {}
    for part in raw.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(None, 1)
            value = value.strip().strip('"')
        else:
            key, value = part, ""
        values = [unquote(item.strip()) for item in value.split(",") if item.strip()]
        attributes.setdefault(key, []).extend(values or [""])
    return {key: tuple(values) for key, values in attributes.items()}


def read_gff3(path: Path) -> tuple[list[Feature], dict[str, str], dict[str, tuple[int, int]]]:
    features: list[Feature] = []
    directives: dict[str, str] = {}
    sequence_regions: dict[str, tuple[int, int]] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##sequence-region "):
                parts = line.split()
                if len(parts) == 4:
                    sequence_regions[parts[1]] = (int(parts[2]), int(parts[3]))
                continue
            if line.startswith("#!"):
                key, _, value = line[2:].partition(" ")
                directives[key] = value.strip()
                continue
            if line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) != 9:
                raise ContractError(f"annotation line {line_number} does not have 9 columns")
            seqid, source, feature_type, start_raw, end_raw, _score, strand, _phase, attrs = columns
            try:
                start_1 = int(start_raw)
                end_1 = int(end_raw)
            except ValueError as exc:
                raise ContractError(f"non-integer coordinate at annotation line {line_number}") from exc
            if start_1 < 1 or end_1 < start_1:
                raise ContractError(f"invalid closed interval at annotation line {line_number}")
            if strand not in {"+", "-", ".", "?"}:
                raise ContractError(f"invalid strand at annotation line {line_number}: {strand}")
            features.append(
                Feature(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start=start_1 - 1,
                    end=end_1,
                    strand=strand,
                    attributes=_parse_attributes(attrs),
                    line_number=line_number,
                )
            )
    if not features:
        raise ContractError("annotation contains no feature records")
    return features, directives, sequence_regions


def feature_type_counts(features: list[Feature]) -> dict[str, int]:
    return dict(sorted(Counter(feature.feature_type for feature in features).items()))

