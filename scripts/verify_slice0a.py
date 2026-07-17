from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import TextIO


TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "snRNA", "snoRNA"}
NCRNA_TYPES = {"ncRNA_gene", "ncRNA", "lnc_RNA", "snRNA", "snoRNA", "miRNA", "SRP_RNA", "RNase_P_RNA", "RNase_MRP_RNA"}
STATISTIC_KEYS = (
    "gene_count",
    "transcript_count",
    "orf_cds_count",
    "ncrna_count",
    "trna_count",
    "rrna_count",
    "partial_gene_count",
    "start_range_gene_count",
    "end_range_gene_count",
    "sequence_count",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def summarize_gff(path: Path) -> dict[str, int]:
    type_counts: Counter[str] = Counter()
    sequence_ids: set[str] = set()
    partial = start_range = end_range = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("##sequence-region "):
                columns = line.split()
                if len(columns) >= 2:
                    sequence_ids.add(columns[1])
                continue
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) != 9:
                raise SystemExit(f"malformed GFF row during independent statistics check: {path}")
            sequence_ids.add(columns[0])
            feature_type = columns[2]
            type_counts[feature_type] += 1
            if feature_type == "gene":
                attributes = {}
                for field in columns[8].split(";"):
                    key, separator, value = field.partition("=")
                    if separator:
                        attributes[key] = value
                partial += attributes.get("partial", "").lower() == "true"
                start_range += "start_range" in attributes
                end_range += "end_range" in attributes
    return {
        "gene_count": type_counts["gene"],
        "transcript_count": sum(type_counts[key] for key in TRANSCRIPT_TYPES),
        "orf_cds_count": type_counts["CDS"],
        "ncrna_count": sum(type_counts[key] for key in NCRNA_TYPES),
        "trna_count": type_counts["tRNA"],
        "rrna_count": type_counts["rRNA"],
        "partial_gene_count": partial,
        "start_range_gene_count": start_range,
        "end_range_gene_count": end_range,
        "sequence_count": len(sequence_ids),
    }


def summarize_kegg(path: Path) -> dict[str, int | None]:
    records = 0
    chromosomes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) < 3:
                raise SystemExit(f"malformed KEGG row during independent statistics check: {path}")
            records += 1
            location = columns[2]
            if ":" in location:
                chromosomes.add(location.split(":", 1)[0].replace("complement(", ""))
    return {
        "gene_count": records,
        "transcript_count": 0,
        "orf_cds_count": records,
        "ncrna_count": 0,
        "trna_count": 0,
        "rrna_count": 0,
        "partial_gene_count": None,
        "start_range_gene_count": None,
        "end_range_gene_count": None,
        "sequence_count": len(chromosomes),
    }


def parse_attributes(value: str) -> dict[str, str]:
    result = {}
    for field in value.split(";"):
        key, separator, item = field.partition("=")
        if separator:
            result[key] = item
    return result


def read_gff_genes(path: Path, sequence_map: dict[str, str] | None = None) -> dict[str, tuple[str, int, int, str]]:
    rows = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) != 9 or columns[2] != "gene":
                continue
            attributes = parse_attributes(columns[8])
            identifier = attributes.get("locus_tag") or attributes.get("gene_id") or attributes.get("ID") or ""
            identifier = identifier.removeprefix("gene:").removeprefix("gene-")
            seqid = (sequence_map or {}).get(columns[0], columns[0])
            rows[identifier] = (seqid, int(columns[3]), int(columns[4]), columns[6])
    return rows


def read_kegg_genes(path: Path, sequence_map: dict[str, str]) -> dict[str, tuple[str, int, int, str]]:
    rows = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            columns = raw.rstrip("\r\n").split("\t")
            if len(columns) < 3:
                continue
            identifier = columns[0].split(":", 1)[-1]
            chromosome, position = columns[2].split(":", 1)
            segments = [(int(start), int(end)) for start, end in re.findall(r"(\d+)\.\.(\d+)", position)]
            if segments:
                rows[identifier] = (
                    sequence_map.get(chromosome, chromosome),
                    min(start for start, _end in segments),
                    max(end for _start, end in segments),
                    "-" if position.startswith("complement(") else "+",
                )
    return rows


def exact_shared_count(left: dict[str, tuple[str, int, int, str]], right: dict[str, tuple[str, int, int, str]]) -> tuple[int, int]:
    shared = set(left) & set(right)
    return len(shared), sum(left[key] == right[key] for key in shared)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run = json.loads((args.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    verified_artifacts = {}
    for name, identity in run["artifacts"].items():
        path = args.run_dir / name
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"artifact identity mismatch: {name}")
        verified_artifacts[name] = identity["sha256"]
    sources = json.loads((args.run_dir / "annotation_sources.json").read_text(encoding="utf-8"))
    qualification = json.loads((args.run_dir / "annotation_qualification.json").read_text(encoding="utf-8"))
    mapping = json.loads((args.run_dir / "annotation_mapping_probe.json").read_text(encoding="utf-8"))
    source_ids = [item["source_id"] for item in sources["sources"]]
    qualification_ids = [item["source_id"] for item in qualification["sources"]]
    if len(source_ids) != len(set(source_ids)) or sorted(source_ids) != sorted(qualification_ids):
        raise SystemExit("source/qualification identity mismatch")
    verified_source_files = 0
    qualification_by_id = {item["source_id"]: item for item in qualification["sources"]}
    source_statistics = {}
    for source in sources["sources"]:
        for spec in source.get("files", {}).values():
            path = args.repo_root / spec["path"]
            if not path.is_file() or path.stat().st_size != spec["size_bytes"] or sha256(path) != spec["sha256"]:
                raise SystemExit(f"source file identity mismatch: {source['source_id']}:{spec['path']}")
            verified_source_files += 1
        statistics = qualification_by_id[source["source_id"]].get("statistics")
        if statistics is None:
            continue
        if source.get("files", {}).get("annotation"):
            recalculated = summarize_gff(args.repo_root / source["files"]["annotation"]["path"])
        elif source.get("annotation_format") == "kegg_list" and source.get("files", {}).get("metadata"):
            recalculated = summarize_kegg(args.repo_root / source["files"]["metadata"]["path"])
        else:
            raise SystemExit(f"no independent statistics route for source: {source['source_id']}")
        expected = {key: statistics.get(key) for key in STATISTIC_KEYS}
        if recalculated != expected:
            raise SystemExit(f"source statistics mismatch: {source['source_id']}")
        source_statistics[source["source_id"]] = recalculated
    source_by_id = {item["source_id"]: item for item in sources["sources"]}
    old_source = source_by_id["strain-b_refseq_gcf000027005_1"]
    ensembl_source = source_by_id["ensembl_fungi_release63_strain-b"]
    kegg_source = source_by_id["kegg_ppa_metadata"]
    old_genes = read_gff_genes(args.repo_root / old_source["files"]["annotation"]["path"])
    ensembl_genes = read_gff_genes(
        args.repo_root / ensembl_source["files"]["annotation"]["path"],
        ensembl_source["sequence_name_map_to_gcf000027005_1"],
    )
    kegg_genes = read_kegg_genes(
        args.repo_root / kegg_source["files"]["metadata"]["path"],
        kegg_source["sequence_name_map_to_gcf000027005_1"],
    )
    ensembl_shared, ensembl_exact = exact_shared_count(old_genes, ensembl_genes)
    kegg_shared, kegg_exact = exact_shared_count(old_genes, kegg_genes)
    mirror_checks = {
        "ensembl_shared": ensembl_shared,
        "ensembl_exact": ensembl_exact,
        "kegg_shared": kegg_shared,
        "kegg_exact": kegg_exact,
    }
    mirror_evidence = qualification["mirror_evidence"]
    if mirror_checks != {
        "ensembl_shared": mirror_evidence["ensembl_vs_old_refseq"]["shared_gene_id_count"],
        "ensembl_exact": mirror_evidence["ensembl_vs_old_refseq"]["exact_boundary_match_count"],
        "kegg_shared": mirror_evidence["kegg_vs_old_refseq"]["shared_gene_id_count"],
        "kegg_exact": mirror_evidence["kegg_vs_old_refseq"]["exact_coordinate_match_count"],
    }:
        raise SystemExit("independent mirror identity mismatch")
    records = mapping["records"]
    summary = mapping["summary"]
    counts = Counter(item["mapping_status"] for item in records)
    for status in ("consistent", "conflict", "unmappable", "uncertain"):
        counts.setdefault(status, 0)
    if dict(sorted(counts.items())) != summary["mapping_status_counts"]:
        raise SystemExit("mapping status summary mismatch")
    order_counts = Counter(item["order_status"] for item in records)
    adjacency_counts = Counter(item["adjacency_status"] for item in records)
    window_order_counts = Counter()
    window_adjacency_counts = Counter()
    for window_id in sorted({item["window_id"] for item in records}):
        window = [item for item in records if item["window_id"] == window_id]
        window_order_counts[window[0]["order_status"]] += 1
        window_adjacency_counts[window[0]["adjacency_status"]] += 1
    if dict(sorted(window_order_counts.items())) != summary["order_status_counts"]:
        raise SystemExit("order status summary mismatch")
    if dict(sorted(window_adjacency_counts.items())) != summary["adjacency_status_counts"]:
        raise SystemExit("adjacency status summary mismatch")
    if set(adjacency_counts) != {"unavailable"}:
        raise SystemExit("mapping probe overclaims target gene adjacency")
    expected_sequences = ["NC_012963.1", "NC_012964.1", "NC_012965.1", "NC_012966.1"]
    if summary["source_nuclear_sequences_covered"] != expected_sequences:
        raise SystemExit("mapping probe does not cover all four old Strain-B chromosomes")
    if ensembl_shared != 5040 or ensembl_exact != 5040:
        raise SystemExit("Ensembl mirror statistics do not prove old-annotation identity")
    if kegg_shared != 5040 or kegg_exact != 5040:
        raise SystemExit("KEGG mirror statistics do not prove old-annotation identity")
    result = {
        "schema_version": 1,
        "verification_type": "slice0a_independent_qualification",
        "status": "passed",
        "run_id": run["run_id"],
        "source_count": len(source_ids),
        "verified_source_file_count": verified_source_files,
        "source_statistics": source_statistics,
        "probe_record_count": len(records),
        "mapping_status_counts": dict(sorted(counts.items())),
        "order_status_counts": dict(sorted(window_order_counts.items())),
        "adjacency_status_counts": dict(sorted(window_adjacency_counts.items())),
        "source_nuclear_sequences_covered": expected_sequences,
        "mirror_checks": mirror_checks,
        "verified_artifacts": verified_artifacts,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
