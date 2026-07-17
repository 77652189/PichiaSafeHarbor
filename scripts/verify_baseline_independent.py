"""Independent Strain-B baseline verifier.

This script intentionally does not import the project package. It reparses FASTA and
top-level GFF gene records, rebuilds boundary unions, checks every published nuclear
intergenic coordinate, verifies chromosome coverage conservation, and records the
first/middle/last interval sampled on each chromosome.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def fasta_lengths(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    name = None
    length = 0
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    result[name] = length
                name = line[1:].split(None, 1)[0]
                length = 0
            elif line:
                length += len(line)
    if name is not None:
        result[name] = length
    return result


def attribute_id(raw: str) -> str:
    for part in raw.split(";"):
        if part.startswith("ID="):
            return part[3:]
    raise ValueError(f"gene record has no ID: {raw}")


def read_genes(path: Path) -> dict[str, list[dict]]:
    result: dict[str, list[dict]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9 or columns[2] != "gene":
                continue
            record = {
                "id": attribute_id(columns[8]),
                "start": int(columns[3]) - 1,
                "end": int(columns[4]),
                "strand": columns[6],
            }
            result.setdefault(columns[0], []).append(record)
    return result


def merge(records: list[dict]) -> list[dict]:
    ordered = sorted(records, key=lambda item: (item["start"], item["end"], item["id"]))
    clusters: list[list[dict]] = []
    for record in ordered:
        if clusters and record["start"] < max(item["end"] for item in clusters[-1]):
            clusters[-1].append(record)
        else:
            clusters.append([record])
    result = []
    for members in clusters:
        start = min(item["start"] for item in members)
        end = max(item["end"] for item in members)
        left = [item for item in members if item["start"] == start]
        right = [item for item in members if item["end"] == end]
        result.append({"start": start, "end": end, "members": members, "left": left, "right": right})
    return result


def orientation(left: str, right: str) -> str:
    if left not in {"+", "-"} or right not in {"+", "-"}:
        return "unknown"
    if left == right:
        return "tandem"
    return "convergent" if (left, right) == ("+", "-") else "divergent"


def strand(records: list[dict]) -> str:
    values = {item["strand"] for item in records if item["strand"] in {"+", "-"}}
    return next(iter(values)) if len(values) == 1 else "."


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--intergenic", type=Path, required=True)
    parser.add_argument("--terminal", type=Path, required=True)
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    run_manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))

    manifest = json.loads(args.reference_manifest.read_text(encoding="utf-8"))["references"]["strain-b"]
    nuclear = sorted(
        seqid for seqid, seq_class in manifest["sequence_classes"].items()
        if seq_class == "nuclear_chromosome"
    )
    lengths = fasta_lengths(args.fasta)
    genes = read_genes(args.gff)
    with args.intergenic.open("r", encoding="utf-8", newline="") as handle:
        published = list(csv.DictReader(handle, delimiter="\t"))
    with args.terminal.open("r", encoding="utf-8", newline="") as handle:
        terminals = list(csv.DictReader(handle, delimiter="\t"))

    expected_all = []
    coverage_checks = {}
    samples = {}
    for seqid in nuclear:
        clusters = merge(genes[seqid])
        expected = []
        for index, (left, right) in enumerate(zip(clusters, clusters[1:]), 1):
            if right["start"] <= left["end"]:
                continue
            left_strand = strand(left["right"])
            right_strand = strand(right["left"])
            expected.append(
                {
                    "seqid": seqid,
                    "start": left["end"],
                    "end": right["start"],
                    "length": right["start"] - left["end"],
                    "left_entity_ids": ",".join(item["id"] for item in left["right"]),
                    "right_entity_ids": ",".join(item["id"] for item in right["left"]),
                    "left_strand": left_strand,
                    "right_strand": right_strand,
                    "orientation": orientation(left_strand, right_strand),
                }
            )
        expected_all.extend(expected)
        actual = [item for item in published if item["seqid"] == seqid]
        comparable_actual = [
            {
                "seqid": item["seqid"],
                "start": int(item["start"]),
                "end": int(item["end"]),
                "length": int(item["length"]),
                "left_entity_ids": item["left_entity_ids"],
                "right_entity_ids": item["right_entity_ids"],
                "left_strand": item["left_strand"],
                "right_strand": item["right_strand"],
                "orientation": item["orientation"],
            }
            for item in actual
        ]
        if comparable_actual != expected:
            raise SystemExit(f"published intergenic records differ from independent calculation: {seqid}")
        terminal_bp = sum(
            int(item["length"]) for item in terminals if item["seqid"] == seqid
        )
        cluster_bp = sum(item["end"] - item["start"] for item in clusters)
        gap_bp = sum(item["length"] for item in expected)
        total = terminal_bp + cluster_bp + gap_bp
        if total != lengths[seqid]:
            raise SystemExit(f"coverage conservation failed for {seqid}: {total} != {lengths[seqid]}")
        coverage_checks[seqid] = {
            "sequence_length": lengths[seqid],
            "boundary_union_bp": cluster_bp,
            "intergenic_bp": gap_bp,
            "terminal_bp": terminal_bp,
            "conserved": True,
        }
        indices = sorted({0, len(expected) // 2, len(expected) - 1})
        samples[seqid] = [expected[index] for index in indices]

    if {item["seqid"] for item in published} != set(nuclear):
        raise SystemExit("published intergenic output contains non-nuclear sequences")
    result = {
        "schema_version": 1,
        "verification_type": "independent_interval_recalculation",
        "run_id": run_manifest["run_id"],
        "status": "passed",
        "method": "independent stdlib FASTA/GFF gene parser; no project-package imports",
        "nuclear_sequences": nuclear,
        "verified_intergenic_region_count": len(expected_all),
        "coverage_checks": coverage_checks,
        "first_middle_last_samples": samples,
        "verified_artifacts": {
            "intergenic_regions.tsv": hashlib.sha256(args.intergenic.read_bytes()).hexdigest(),
            "terminal_regions.tsv": hashlib.sha256(args.terminal.read_bytes()).hexdigest(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
