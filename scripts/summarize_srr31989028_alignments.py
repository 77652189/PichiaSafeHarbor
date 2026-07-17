from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
NUCLEAR = {"CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"}
CANONICAL_UNSTRANDED = {"GT-AG", "GC-AG", "AT-AC", "CT-AC", "CT-GC", "GT-AT"}


def load_fasta(path: Path) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            result[current] = []
        elif current:
            result[current].append(line.strip())
    return {key: "".join(value).upper() for key, value in result.items()}


def tag_int(fields: list[str], name: str) -> int | None:
    prefix = name + ":i:"
    for value in fields[11:]:
        if value.startswith(prefix):
            return int(value[len(prefix):])
    return None


def junctions(rname: str, pos1: int, cigar: str) -> list[tuple[str, int, int]]:
    ref = pos1 - 1
    result = []
    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op == "N":
            result.append((rname, ref, ref + length))
            ref += length
        elif op in {"M", "D", "=", "X"}:
            ref += length
    return result


def process_group(records: list[list[str]], counts: Counter, junction_counts: dict) -> None:
    if not records:
        return
    counts["total_query_pairs"] += 1
    primary = [row for row in records if not (int(row[1]) & (0x100 | 0x800))]
    r1 = next((row for row in primary if int(row[1]) & 0x40), None)
    r2 = next((row for row in primary if int(row[1]) & 0x80), None)
    if r1 is None or r2 is None:
        counts["unknown_pair_structure"] += 1
    else:
        mapped = [not (int(row[1]) & 0x4) for row in (r1, r2)]
        if mapped == [False, False]:
            counts["pair_unmapped"] += 1
        elif mapped[0] != mapped[1]:
            counts["pair_mixed"] += 1
        else:
            counts["pair_both_mapped"] += 1
            if (int(r1[1]) & 0x2) and (int(r2[1]) & 0x2):
                counts["pair_proper"] += 1
            else:
                counts["pair_discordant"] += 1
            nh = [tag_int(row, "NH") for row in (r1, r2)]
            has_secondary = any(int(row[1]) & 0x100 for row in records)
            if nh == [1, 1] and not has_secondary:
                counts["pair_unique"] += 1
            elif has_secondary or any(value is not None and value > 1 for value in nh):
                counts["pair_multi"] += 1
            else:
                counts["pair_multiplicity_unknown"] += 1
        if any(row[5] != "*" and any(op == "S" and int(length) >= 20 for length, op in CIGAR_RE.findall(row[5])) for row in (r1, r2)):
            counts["pair_softclip_ge20"] += 1
        if any(row[5] != "*" and any(op in {"I", "D"} and int(length) >= 20 for length, op in CIGAR_RE.findall(row[5])) for row in (r1, r2)):
            counts["pair_indel_ge20"] += 1
        mapped_refs = {row[2] for row in (r1, r2) if row[2] != "*"}
        if mapped_refs and not mapped_refs.issubset(NUCLEAR):
            counts["pair_non_nuclear_reference"] += 1
    if any(int(row[1]) & 0x800 for row in records):
        counts["pair_with_supplementary"] += 1
    unique_seen = set(); multi_seen = set()
    for row in records:
        flag = int(row[1])
        if flag & (0x4 | 0x800) or row[2] == "*" or row[5] == "*":
            continue
        nh = tag_int(row, "NH")
        target = unique_seen if nh == 1 and not (flag & (0x100 | 0x800)) else multi_seen
        target.update(junctions(row[2], int(row[3]), row[5]))
    for key in unique_seen:
        junction_counts[key]["unique"] += 1
    for key in multi_seen:
        junction_counts[key]["multi"] += 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    parser.add_argument("--junctions-tsv", type=Path, required=True)
    args = parser.parse_args()
    sequences = load_fasta(args.reference_fasta)
    counts: Counter = Counter()
    junction_counts = defaultdict(lambda: {"unique": 0, "multi": 0})
    current = None; group = []
    for line in sys.stdin:
        if not line or line.startswith("@"):
            continue
        row = line.rstrip("\n").split("\t")
        counts["alignment_records"] += 1
        flag = int(row[1])
        if flag & 0x100: counts["secondary_alignment_records"] += 1
        if flag & 0x800: counts["supplementary_alignment_records"] += 1
        if current is not None and row[0] != current:
            process_group(group, counts, junction_counts)
            group = []
        current = row[0]; group.append(row)
    process_group(group, counts, junction_counts)
    total = counts["total_query_pairs"]
    classified = counts["pair_unmapped"] + counts["pair_mixed"] + counts["pair_both_mapped"] + counts["unknown_pair_structure"]
    summary = {
        "schema_version": 1,
        "classification_unit": "query_pair grouped after samtools collate; unique/multi applies only to pairs with both mates mapped",
        "strand_support": "unavailable",
        "counts": dict(sorted(counts.items())),
        "classification_complete": classified == total,
        "junction_summary": {
            "distinct_junctions": len(junction_counts),
            "with_unique_fragment_support": sum(value["unique"] > 0 for value in junction_counts.values()),
            "with_multi_fragment_support": sum(value["multi"] > 0 for value in junction_counts.values()),
        },
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with args.junctions_tsv.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("reference\tintron_start_0based\tintron_end_0based\tmotif_unstranded\tcanonical_unstranded\tunique_fragment_support\tmulti_fragment_support\tstrand_support\n")
        for (reference, start, end), support in sorted(junction_counts.items()):
            seq = sequences.get(reference, "")
            motif = f"{seq[start:start+2]}-{seq[end-2:end]}" if 0 <= start < end <= len(seq) else "unavailable"
            handle.write(f"{reference}\t{start}\t{end}\t{motif}\t{str(motif in CANONICAL_UNSTRANDED).lower()}\t{support['unique']}\t{support['multi']}\tunavailable\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
