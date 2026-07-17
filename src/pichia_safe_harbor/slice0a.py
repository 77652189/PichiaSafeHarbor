from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

from .annotation_qualification import (
    compare_gff_gene_mirror,
    compare_kegg_to_gff,
    load_annotation_source_manifest,
    summarize_sources,
)
from .errors import ContractError
from .io_utils import sha256_file, write_json, write_tsv
from .mapping_probe import run_mapping_probe
from .pipeline import _implementation_hash


def _source(manifest: dict[str, Any], source_id: str) -> dict[str, Any]:
    return next(source for source in manifest["sources"] if source["source_id"] == source_id)


def _run_id(manifest_hash: str, evidence_hash: str, implementation_hash: str) -> str:
    material = json.dumps(
        {
            "source_manifest_sha256": manifest_hash,
            "acquisition_evidence_sha256": evidence_hash,
            "implementation_sha256": implementation_hash,
            "mapping_method": "exact-101bp-boundary-anchors-v1",
            "primary_coordinate_space": "GCA_001746955.1",
        },
        sort_keys=True,
    ).encode("utf-8")
    return "slice0a-" + hashlib.sha256(material).hexdigest()[:16]


def run_slice0a(source_manifest_path: Path, output_dir: Path, repo_root: Path) -> dict[str, Any]:
    if output_dir.exists():
        raise ContractError(f"output directory already exists: {output_dir}")
    manifest = load_annotation_source_manifest(source_manifest_path)
    evidence_path = repo_root / manifest["acquisition_evidence"]
    if not evidence_path.is_file():
        raise ContractError(f"acquisition evidence is missing: {evidence_path}")
    qualifications = summarize_sources(manifest, repo_root)

    old = _source(manifest, "strain-b_refseq_gcf000027005_1")
    current = _source(manifest, "strain-b_insdc_gca001746955_1")
    ensembl = _source(manifest, "ensembl_fungi_release63_strain-b")
    kegg = _source(manifest, "kegg_ppa_metadata")
    old_gff = repo_root / old["files"]["annotation"]["path"]
    mirror_evidence = {
        "ensembl_vs_old_refseq": compare_gff_gene_mirror(
            old_gff,
            repo_root / ensembl["files"]["annotation"]["path"],
            ensembl["sequence_name_map_to_gcf000027005_1"],
        ),
        "kegg_vs_old_refseq": compare_kegg_to_gff(
            repo_root / kegg["files"]["metadata"]["path"],
            old_gff,
            kegg["sequence_name_map_to_gcf000027005_1"],
        ),
    }
    for item in qualifications:
        if item["source_id"] == "ensembl_fungi_release63_strain-b":
            item["mirror_evidence"] = mirror_evidence["ensembl_vs_old_refseq"]
        elif item["source_id"] == "kegg_ppa_metadata":
            item["mirror_evidence"] = mirror_evidence["kegg_vs_old_refseq"]

    mapping_records, mapping_summary = run_mapping_probe(
        repo_root / old["files"]["fasta"]["path"],
        old_gff,
        repo_root / current["files"]["fasta"]["path"],
        {"CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"},
    )

    source_manifest_hash = sha256_file(source_manifest_path)
    evidence_hash = sha256_file(evidence_path)
    implementation_hash = _implementation_hash()
    run_id = _run_id(source_manifest_hash, evidence_hash, implementation_hash)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        acquisition_evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
        source_output = {
            "schema_version": 1,
            "manifest_version": manifest["manifest_version"],
            "primary_coordinate_space": manifest["primary_coordinate_space"],
            "source_manifest": {"path": str(source_manifest_path), "size_bytes": source_manifest_path.stat().st_size, "sha256": source_manifest_hash},
            "acquisition_evidence": {
                "path": manifest["acquisition_evidence"],
                "size_bytes": evidence_path.stat().st_size,
                "sha256": evidence_hash,
                "checked_on": acquisition_evidence["checked_on"],
                "timezone": acquisition_evidence["timezone"],
                "probes": acquisition_evidence["probes"],
            },
            "sources": manifest["sources"],
        }
        write_json(temp_dir / "annotation_sources.json", source_output)
        write_json(temp_dir / "annotation_qualification.json", {"schema_version": 1, "sources": qualifications, "mirror_evidence": mirror_evidence})
        qualification_rows = []
        for item in qualifications:
            stats = item.get("statistics") or {}
            qualification_rows.append(
                {
                    "source_id": item["source_id"],
                    "availability": item["availability"],
                    "source_version": item["source_version"],
                    "source_url": item["source_url"],
                    "strain": item["strain"],
                    "assembly_accession": item["assembly_accession"],
                    "coordinate_space": item["coordinate_space"],
                    "qualification_status": item["qualification_status"],
                    "gene_count": stats.get("gene_count"),
                    "transcript_count": stats.get("transcript_count"),
                    "orf_cds_count": stats.get("orf_cds_count"),
                    "ncrna_count": stats.get("ncrna_count"),
                    "trna_count": stats.get("trna_count"),
                    "rrna_count": stats.get("rrna_count"),
                    "partial_gene_count": stats.get("partial_gene_count"),
                    "partial_gene_fraction": stats.get("partial_gene_fraction"),
                    "sequence_count": stats.get("sequence_count"),
                    "independence_group": item["upstream"]["independence_group"],
                    "upstream_relationship": item["upstream"]["relationship"],
                    "qualification_reasons": item["qualification_reasons"],
                    "license_status": item["license"]["status"],
                    "license_url": item["license"]["url"],
                }
            )
        write_tsv(temp_dir / "annotation_qualification.tsv", qualification_rows, list(qualification_rows[0]))
        mapping_output = {"schema_version": 1, "source_assembly": old["assembly_accession"], "target_assembly": "GCA_001746955.1", "records": mapping_records, "summary": mapping_summary}
        write_json(temp_dir / "annotation_mapping_probe.json", mapping_output)
        mapping_fields = sorted({key for record in mapping_records for key in record})
        write_tsv(temp_dir / "annotation_mapping_probe.tsv", mapping_records, mapping_fields)

        statuses = {item["source_id"]: item["qualification_status"] for item in qualifications}
        report = _render_report(qualifications, mapping_summary, mirror_evidence)
        (temp_dir / "annotation_qualification_report.md").write_text(report, encoding="utf-8")
        recommendation = _render_recommendation(statuses, mapping_summary)
        (temp_dir / "source_selection_recommendation.md").write_text(recommendation, encoding="utf-8")
        artifacts = [
            "annotation_sources.json",
            "annotation_qualification.json",
            "annotation_qualification.tsv",
            "annotation_mapping_probe.json",
            "annotation_mapping_probe.tsv",
            "annotation_qualification_report.md",
            "source_selection_recommendation.md",
        ]
        run_manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "slice": "Slice 0A annotation source qualification",
            "execution_status": "complete",
            "verification_status": "not_run",
            "scientific_acceptance_status": "blocked",
            "scientific_acceptance_blockers": ["no qualified direct annotation source", "new source-selection ADR required"],
            "primary_coordinate_space": "GCA_001746955.1",
            "implementation_sha256": implementation_hash,
            "source_manifest_sha256": source_manifest_hash,
            "acquisition_evidence_sha256": evidence_hash,
            "mapping_configuration": {"anchor_length": 101, "selection": "first/middle/last consecutive-three-gene windows on all four old Strain-B chromosomes"},
            "artifacts": {name: {"size_bytes": (temp_dir / name).stat().st_size, "sha256": sha256_file(temp_dir / name)} for name in artifacts},
        }
        write_json(temp_dir / "run_manifest.json", run_manifest)
        os.replace(temp_dir, output_dir)
        return run_manifest
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


def create_slice0a_acceptance(run_dir: Path, independent_path: Path, test_evidence_path: Path) -> dict[str, Any]:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    independent = json.loads(independent_path.read_text(encoding="utf-8"))
    tests = json.loads(test_evidence_path.read_text(encoding="utf-8"))
    expected_artifacts = {
        "annotation_sources.json",
        "annotation_qualification.json",
        "annotation_qualification.tsv",
        "annotation_mapping_probe.json",
        "annotation_mapping_probe.tsv",
        "annotation_qualification_report.md",
        "source_selection_recommendation.md",
    }
    if set(run.get("artifacts", {})) != expected_artifacts:
        raise ContractError("Slice 0A artifact set does not satisfy the completion contract")
    if run.get("primary_coordinate_space") != "GCA_001746955.1" or run.get("scientific_acceptance_status") != "blocked":
        raise ContractError("Slice 0A run scope or scientific status mismatch")
    if independent.get("schema_version") != 1 or independent.get("verification_type") != "slice0a_independent_qualification" or independent.get("status") != "passed" or independent.get("run_id") != run["run_id"]:
        raise ContractError("independent Slice 0A evidence did not pass its identity contract")
    if tests.get("schema_version") != 1 or tests.get("evidence_type") != "automated_tests" or tests.get("status") != "passed" or tests.get("implementation_sha256") != run["implementation_sha256"] or tests.get("passed_count", 0) <= 0 or tests.get("command") != ["python", "-m", "pytest", "-q"]:
        raise ContractError("automated test evidence did not pass its identity contract")
    for name, identity in run["artifacts"].items():
        path = run_dir / name
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256_file(path) != identity["sha256"]:
            raise ContractError(f"Slice 0A artifact identity mismatch: {name}")
    if independent.get("verified_artifacts") != {
        name: identity["sha256"] for name, identity in run["artifacts"].items()
    }:
        raise ContractError("independent Slice 0A artifact identity summary mismatch")
    sources = json.loads((run_dir / "annotation_sources.json").read_text(encoding="utf-8"))
    qualification = json.loads((run_dir / "annotation_qualification.json").read_text(encoding="utf-8"))
    mapping = json.loads((run_dir / "annotation_mapping_probe.json").read_text(encoding="utf-8"))
    required_source_ids = {
        "strain-b_insdc_gca001746955_1",
        "strain-b_refseq_gcf000027005_1",
        "strain-c_genbank_gca000223565_1",
        "ensembl_fungi_release63_strain-b",
        "kegg_ppa_metadata",
        "pichiagenome_2016_curation",
        "biocyc_picpa",
    }
    source_by_id = {item["source_id"]: item for item in sources["sources"]}
    if set(source_by_id) != required_source_ids or sources.get("primary_coordinate_space") != "GCA_001746955.1":
        raise ContractError("Slice 0A required source inventory mismatch")
    for source in source_by_id.values():
        if source["availability"] == "unavailable":
            continue
        if not source.get("source_version") or not source.get("source_url") or not source.get("coordinate_space"):
            raise ContractError(f"available source identity is incomplete: {source['source_id']}")
        license_record = source.get("license") or {}
        if not license_record.get("status") or not license_record.get("url"):
            raise ContractError(f"available source license is incomplete: {source['source_id']}")
        if not source.get("files"):
            raise ContractError(f"available source has no hashed files: {source['source_id']}")
        for role, identity in source["files"].items():
            if not identity.get("path") or not identity.get("sha256") or not isinstance(identity.get("size_bytes"), int):
                raise ContractError(f"available source file identity is incomplete: {source['source_id']}:{role}")
    pichiagenome = source_by_id["pichiagenome_2016_curation"]
    if pichiagenome["availability"] != "unavailable" or not {"RNA-seq", "proteomics"}.issubset(pichiagenome.get("experimental_support", [])):
        raise ContractError("PichiaGenome qualification record is incomplete")
    probes = sources.get("acquisition_evidence", {}).get("probes", [])
    pichia_probes = [item for item in probes if item.get("probe_id", "").startswith("pichiagenome_")]
    if not pichia_probes or not any(item.get("availability") == "unavailable" for item in pichia_probes):
        raise ContractError("PichiaGenome reproducible acquisition failure evidence is missing")
    if any(not item.get("target") or not item.get("method") or not item.get("result") for item in pichia_probes):
        raise ContractError("PichiaGenome acquisition failure evidence is incomplete")
    old = source_by_id["strain-b_refseq_gcf000027005_1"]
    for mirror_id in ("ensembl_fungi_release63_strain-b", "kegg_ppa_metadata"):
        mirror = source_by_id[mirror_id]
        if mirror["upstream"]["independence_group"] != old["upstream"]["independence_group"] or "mirror" not in mirror["upstream"]["relationship"]:
            raise ContractError(f"mirror upstream identity is not explicit: {mirror_id}")
    expected_qualification = {
        "strain-b_insdc_gca001746955_1": "unqualified",
        "strain-b_refseq_gcf000027005_1": "unqualified",
        "strain-c_genbank_gca000223565_1": "unqualified",
        "ensembl_fungi_release63_strain-b": "unqualified",
        "kegg_ppa_metadata": "unqualified",
        "pichiagenome_2016_curation": "unavailable",
        "biocyc_picpa": "unavailable",
    }
    actual_qualification = {item["source_id"]: item["qualification_status"] for item in qualification["sources"]}
    if actual_qualification != expected_qualification:
        raise ContractError("Slice 0A qualification status contract mismatch")
    expected_sequences = ["NC_012963.1", "NC_012964.1", "NC_012965.1", "NC_012966.1"]
    required_mapping_statuses = {"consistent", "conflict", "unmappable", "uncertain"}
    mapping_counts = mapping["summary"]["mapping_status_counts"]
    if set(mapping_counts) != required_mapping_statuses or any(mapping_counts[status] <= 0 for status in required_mapping_statuses):
        raise ContractError("Slice 0A mapping probe does not contain all required classifications")
    if len(mapping["records"]) != 36 or mapping["summary"]["source_nuclear_sequences_covered"] != expected_sequences:
        raise ContractError("Slice 0A mapping probe coverage mismatch")
    for seqid in expected_sequences:
        positions = {item["window_position"] for item in mapping["records"] if item.get("source_seqid") == seqid}
        if positions != {"first", "middle", "last"}:
            raise ContractError(f"Slice 0A mapping windows are incomplete: {seqid}")
    if mapping["summary"].get("adjacency_status_counts") != {"unavailable": 12}:
        raise ContractError("Slice 0A mapping probe overclaims target gene adjacency")
    if independent.get("source_count") != len(sources["sources"]):
        raise ContractError("independent Slice 0A source count mismatch")
    statistic_keys = (
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
    expected_statistics = {
        item["source_id"]: {key: item["statistics"].get(key) for key in statistic_keys}
        for item in qualification["sources"]
        if item.get("statistics") is not None
    }
    if independent.get("source_statistics") != expected_statistics:
        raise ContractError("independent Slice 0A source statistics mismatch")
    if independent.get("probe_record_count") != len(mapping["records"]):
        raise ContractError("independent Slice 0A probe count mismatch")
    if independent.get("mapping_status_counts") != mapping["summary"]["mapping_status_counts"]:
        raise ContractError("independent Slice 0A mapping summary mismatch")
    if independent.get("order_status_counts") != mapping["summary"]["order_status_counts"]:
        raise ContractError("independent Slice 0A order summary mismatch")
    if independent.get("adjacency_status_counts") != mapping["summary"]["adjacency_status_counts"]:
        raise ContractError("independent Slice 0A adjacency summary mismatch")
    if independent.get("source_nuclear_sequences_covered") != mapping["summary"]["source_nuclear_sequences_covered"]:
        raise ContractError("independent Slice 0A chromosome coverage mismatch")
    expected_mirror_checks = {
        "ensembl_shared": qualification["mirror_evidence"]["ensembl_vs_old_refseq"]["shared_gene_id_count"],
        "ensembl_exact": qualification["mirror_evidence"]["ensembl_vs_old_refseq"]["exact_boundary_match_count"],
        "kegg_shared": qualification["mirror_evidence"]["kegg_vs_old_refseq"]["shared_gene_id_count"],
        "kegg_exact": qualification["mirror_evidence"]["kegg_vs_old_refseq"]["exact_coordinate_match_count"],
    }
    if independent.get("mirror_checks") != expected_mirror_checks or set(expected_mirror_checks.values()) != {5040}:
        raise ContractError("independent Slice 0A mirror identity mismatch")
    report = (run_dir / "annotation_qualification_report.md").read_text(encoding="utf-8")
    recommendation = (run_dir / "source_selection_recommendation.md").read_text(encoding="utf-8")
    if "28 high-confidence intervals are not accepted" not in report or "No thresholds, candidate windows, or formal candidates are generated" not in report:
        raise ContractError("Slice 0A report scope guard mismatch")
    if "do not start it without a new ADR" not in recommendation or "Scientific status: blocked" not in recommendation:
        raise ContractError("Slice 0A recommendation stop line mismatch")
    result = {
        "schema_version": 1,
        "acceptance_version": "slice0a-adr0004-v1",
        "run_id": run["run_id"],
        "execution_status": "complete",
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
        "scientific_acceptance_blockers": ["no qualified direct boundary source", "source-selection ADR pending"],
        "run_manifest": {"size_bytes": (run_dir / "run_manifest.json").stat().st_size, "sha256": sha256_file(run_dir / "run_manifest.json")},
        "artifacts": run["artifacts"],
        "verification_evidence": {
            "independent": {"size_bytes": independent_path.stat().st_size, "sha256": sha256_file(independent_path)},
            "automated_tests": {"size_bytes": test_evidence_path.stat().st_size, "sha256": sha256_file(test_evidence_path), "passed_count": tests["passed_count"]},
        },
        "next_authoritative_action": "create a new ADR selecting a boundary source, consensus track, or RNA-seq-supported reannotation route",
    }
    write_json(run_dir / "acceptance_manifest.json", result)
    return result


def _render_report(qualifications: list[dict[str, Any]], mapping: dict[str, Any], mirrors: dict[str, Any]) -> str:
    lines = ["# Slice 0A annotation source qualification report", "", "## Scope guard", "", "- Primary coordinate assembly remains Strain-B GCA_001746955.1.", "- The 28 high-confidence intervals are not accepted as threshold evidence.", "- No thresholds, candidate windows, or formal candidates are generated.", "", "## Source qualification", "", "| Source | Availability | Assembly | Genes | Partial fraction | Independence | Qualification |", "| --- | --- | --- | ---: | ---: | --- | --- |"]
    for item in qualifications:
        stats = item.get("statistics") or {}
        lines.append(f"| {item['source_id']} | {item['availability']} | {item['assembly_accession']} | {stats.get('gene_count', 'unavailable')} | {stats.get('partial_gene_fraction', 'unavailable')} | {item['upstream']['relationship']} | {item['qualification_status']} |")
    lines.extend(["", "## Mirror identity", "", f"- Ensembl/old RefSeq shared IDs: {mirrors['ensembl_vs_old_refseq']['shared_gene_id_count']}; exact chromosome-and-boundary matches: {mirrors['ensembl_vs_old_refseq']['exact_boundary_match_count']}.", f"- KEGG/old RefSeq shared IDs: {mirrors['kegg_vs_old_refseq']['shared_gene_id_count']}; exact chromosome-and-coordinate matches: {mirrors['kegg_vs_old_refseq']['exact_coordinate_match_count']}.", "", "## Four-chromosome mapping probe", "", f"- Records: {mapping['probe_record_count']}", f"- Source chromosomes covered: {mapping['source_nuclear_sequence_count']}", f"- Mapping classes: {json.dumps(mapping['mapping_status_counts'], sort_keys=True)}", f"- Collinear-order classes: {json.dumps(mapping['order_status_counts'], sort_keys=True)}", f"- Gene adjacency classes: {json.dumps(mapping['adjacency_status_counts'], sort_keys=True)}; exact anchors do not qualify target gene adjacency.", "", "## Qualification conclusion", "", "No available source is qualified as a direct scientific boundary source for GCA_001746955.1. PichiaGenome remains the highest-value experimentally supported candidate but its annotation file and exact coordinate identity are unavailable. Scientific status remains blocked pending a new source-selection ADR."])
    return "\n".join(lines) + "\n"


def _render_recommendation(statuses: dict[str, str], mapping: dict[str, Any]) -> str:
    return "\n".join([
        "# Annotation source selection recommendation (pending ADR)",
        "",
        "Recommended route: plan a versioned consensus-boundary or RNA-seq-supported reannotation track on the locked GCA_001746955.1 coordinate space, but do not start it without a new ADR.",
        "",
        "Rationale:",
        "",
        "- Current Strain-B INSDC and old Strain-B RefSeq annotations have systematic partial-boundary flags.",
        "- Ensembl and KEGG trace to the old ASM2700v1 annotation and are not independent boundary evidence.",
        "- Strain-C is a different strain and its current GenBank annotation has the same partial-boundary problem.",
        "- PichiaGenome has RNA-seq/proteomics support but its persistent annotation file could not be recovered.",
        f"- The small mapping probe covers {mapping['source_nuclear_sequence_count']} old Strain-B chromosomes and demonstrates that explicit mapping is feasible for a subset, but it does not authorize whole-genome source selection.",
        "",
        "Alternative if PichiaGenome data are recovered: qualify and map that curated Strain-C track before deciding whether to build a consensus boundary track.",
        "",
        "Scientific status: blocked. Thresholds and formal candidates remain unavailable.",
    ]) + "\n"
