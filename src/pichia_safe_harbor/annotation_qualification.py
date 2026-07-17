from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file
from .parsers import read_gff3

TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "snRNA", "snoRNA"}
NCRNA_TYPES = {"ncRNA_gene", "ncRNA", "lnc_RNA", "snRNA", "snoRNA", "miRNA", "SRP_RNA", "RNase_P_RNA", "RNase_MRP_RNA"}


def load_annotation_source_manifest(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema_version") != 1 or not value.get("sources"):
        raise ContractError("unsupported or empty annotation source manifest")
    return value


def verify_source_files(source: dict[str, Any], repo_root: Path) -> dict[str, dict[str, Any]]:
    verified: dict[str, dict[str, Any]] = {}
    for role, spec in sorted(source.get("files", {}).items()):
        path = repo_root / spec["path"]
        if not path.is_file():
            raise ContractError(f"source file is missing: {source['source_id']}:{role}:{path}")
        actual = {"path": spec["path"], "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
        if actual["size_bytes"] != spec["size_bytes"] or actual["sha256"] != spec["sha256"]:
            raise ContractError(f"source file identity mismatch: {source['source_id']}:{role}")
        verified[role] = actual
    if source["availability"] == "available" and not verified:
        raise ContractError(f"available source has no verifiable files: {source['source_id']}")
    return verified


def _attribute_true(feature, key: str) -> bool:
    return any(value.lower() == "true" for value in feature.attributes.get(key, ()))


def summarize_gff(path: Path) -> dict[str, Any]:
    features, directives, sequence_regions = read_gff3(path)
    type_counts = Counter(feature.feature_type for feature in features)
    genes = [feature for feature in features if feature.feature_type == "gene"]
    partial = [feature for feature in genes if _attribute_true(feature, "partial")]
    start_range = [feature for feature in genes if feature.attributes.get("start_range")]
    end_range = [feature for feature in genes if feature.attributes.get("end_range")]
    return {
        "feature_type_counts": dict(sorted(type_counts.items())),
        "gene_count": len(genes),
        "transcript_count": sum(type_counts.get(key, 0) for key in TRANSCRIPT_TYPES),
        "orf_cds_count": type_counts.get("CDS", 0),
        "ncrna_count": sum(type_counts.get(key, 0) for key in NCRNA_TYPES),
        "trna_count": type_counts.get("tRNA", 0),
        "rrna_count": type_counts.get("rRNA", 0),
        "partial_gene_count": len(partial),
        "start_range_gene_count": len(start_range),
        "end_range_gene_count": len(end_range),
        "partial_gene_fraction": round(len(partial) / len(genes), 6) if genes else None,
        "sequence_count": len(sequence_regions or {feature.seqid for feature in features}),
        "sequence_ids": sorted(sequence_regions or {feature.seqid for feature in features}),
        "directives": directives,
        "missing_feature_types": [
            label for label, count in (("ncRNA", sum(type_counts.get(key, 0) for key in NCRNA_TYPES)), ("tRNA", type_counts.get("tRNA", 0)), ("rRNA", type_counts.get("rRNA", 0))) if count == 0
        ],
    }


def summarize_kegg_list(path: Path) -> dict[str, Any]:
    records = [line.rstrip("\r\n").split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    chromosomes = set()
    for columns in records:
        if len(columns) >= 3:
            match = re.match(r"([^:]+):", columns[2])
            if match:
                chromosomes.add(match.group(1).replace("complement(", ""))
    return {
        "feature_type_counts": {"CDS": len(records)},
        "gene_count": len(records),
        "transcript_count": 0,
        "orf_cds_count": len(records),
        "ncrna_count": 0,
        "trna_count": 0,
        "rrna_count": 0,
        "partial_gene_count": None,
        "start_range_gene_count": None,
        "end_range_gene_count": None,
        "partial_gene_fraction": None,
        "sequence_count": len(chromosomes),
        "sequence_ids": sorted(chromosomes),
        "directives": {},
        "missing_feature_types": ["ncRNA", "tRNA", "rRNA"],
    }


def qualification_status(source: dict[str, Any], stats: dict[str, Any] | None) -> tuple[str, list[str]]:
    if source["availability"] == "unavailable":
        return "unavailable", ["source file could not be obtained or verified"]
    relationship = source["upstream"]["relationship"]
    reasons: list[str] = []
    if "mirror" in relationship or source["expected_use"] == "mirror_detection_only":
        reasons.append("same upstream annotation or coordinate mirror; not independent evidence")
    if source["strain"] != "Strain-B":
        reasons.append("different strain; cannot directly define Strain-B primary-coordinate boundaries")
    if stats and stats.get("partial_gene_fraction") is not None and stats["partial_gene_fraction"] > 0.5:
        reasons.append("systematic partial gene boundaries exceed 50 percent")
    if source["coordinate_space"] != "GCA_001746955.1":
        reasons.append("requires explicit mapping to the locked primary assembly")
    if source["availability"] == "metadata_only":
        reasons.append("metadata does not provide an independently qualified boundary file")
    if reasons:
        return "unqualified", reasons
    return "qualified", ["same primary coordinate space with sufficiently explicit boundaries"]


def summarize_sources(manifest: dict[str, Any], repo_root: Path) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for source in manifest["sources"]:
        verified = verify_source_files(source, repo_root)
        stats: dict[str, Any] | None = None
        annotation = source.get("files", {}).get("annotation")
        metadata = source.get("files", {}).get("metadata")
        if annotation:
            stats = summarize_gff(repo_root / annotation["path"])
        elif source["annotation_format"] == "kegg_list" and metadata:
            stats = summarize_kegg_list(repo_root / metadata["path"])
        status, reasons = qualification_status(source, stats)
        results.append(
            {
                "source_id": source["source_id"],
                "display_name": source["display_name"],
                "source_version": source["source_version"],
                "source_url": source["source_url"],
                "availability": source["availability"],
                "strain": source["strain"],
                "assembly_accession": source["assembly_accession"],
                "coordinate_space": source["coordinate_space"],
                "license": source["license"],
                "upstream": source["upstream"],
                "experimental_support": source["experimental_support"],
                "boundary_uncertainty_semantics": source["boundary_uncertainty_semantics"],
                "verified_files": verified,
                "statistics": stats,
                "qualification_status": status,
                "qualification_reasons": reasons,
            }
        )
    return results


def compare_gff_gene_mirror(
    left_path: Path,
    right_path: Path,
    right_sequence_map: dict[str, str] | None = None,
) -> dict[str, Any]:
    left, _, _ = read_gff3(left_path)
    right, _, _ = read_gff3(right_path)

    def canonical(features, sequence_map: dict[str, str] | None = None):
        rows = set()
        for feature in features:
            if feature.feature_type != "gene":
                continue
            identifier = (
                feature.attributes.get("locus_tag")
                or feature.attributes.get("gene_id")
                or feature.attributes.get("ID")
                or ("",)
            )[0].removeprefix("gene:").removeprefix("gene-")
            seqid = (sequence_map or {}).get(feature.seqid, feature.seqid)
            rows.add((identifier, seqid, feature.start, feature.end, feature.strand))
        return rows

    left_rows = canonical(left)
    right_rows = canonical(right, right_sequence_map)
    common_ids = {row[0] for row in left_rows} & {row[0] for row in right_rows}
    left_by_id = {row[0]: row[1:] for row in left_rows}
    right_by_id = {row[0]: row[1:] for row in right_rows}
    exact_boundaries = sum(left_by_id[key] == right_by_id[key] for key in common_ids)
    return {
        "left_gene_count": len(left_rows),
        "right_gene_count": len(right_rows),
        "shared_gene_id_count": len(common_ids),
        "exact_boundary_match_count": exact_boundaries,
        "exact_boundary_match_fraction": round(exact_boundaries / len(common_ids), 6) if common_ids else 0.0,
    }


def compare_kegg_to_gff(
    kegg_path: Path,
    gff_path: Path,
    kegg_sequence_map: dict[str, str] | None = None,
) -> dict[str, Any]:
    features, _, _ = read_gff3(gff_path)
    gff_rows: dict[str, tuple[str, int, int, str]] = {}
    for feature in features:
        if feature.feature_type != "gene":
            continue
        values = feature.attributes.get("locus_tag") or feature.attributes.get("Name")
        if values:
            gff_rows[values[0]] = (feature.seqid, feature.start + 1, feature.end, feature.strand)
    kegg_rows: dict[str, tuple[str, int, int, str]] = {}
    for line in kegg_path.read_text(encoding="utf-8").splitlines():
        columns = line.split("\t")
        if len(columns) < 3:
            continue
        identifier = columns[0].split(":", 1)[-1]
        position = columns[2].split(":", 1)[-1]
        chromosome = columns[2].split(":", 1)[0].replace("complement(", "")
        seqid = (kegg_sequence_map or {}).get(chromosome, chromosome)
        strand = "-" if position.startswith("complement(") else "+"
        segments = [(int(start), int(end)) for start, end in re.findall(r"(\d+)\.\.(\d+)", position)]
        if segments:
            kegg_rows[identifier] = (
                seqid,
                min(start for start, _end in segments),
                max(end for _start, end in segments),
                strand,
            )
    common = sorted(set(kegg_rows) & set(gff_rows))
    exact = sum(kegg_rows[key] == gff_rows[key] for key in common)
    return {
        "kegg_record_count": len(kegg_rows),
        "gff_gene_count": len(gff_rows),
        "shared_gene_id_count": len(common),
        "exact_coordinate_match_count": exact,
        "exact_coordinate_match_fraction": round(exact / len(common), 6) if common else 0.0,
    }
