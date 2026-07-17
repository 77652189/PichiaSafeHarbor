from __future__ import annotations

import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
REQUIRED_CLASSES = ("success", "conflict", "multi_mapping", "unmappable", "endpoint_limited")
SEQUENCE_MAP = {"chr1": "CP014715.1", "chr2": "CP014716.1", "chr3": "CP014717.1", "chr4": "CP014718.1"}
WINDOWS = ("start", "middle", "end")
PROBE_MODE = "submitted_alignment_coordinate_probe"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    sequences = {}
    name = None
    parts = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts).upper()
                name = line[1:].split(None, 1)[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        sequences[name] = "".join(parts).upper()
    return sequences


def recalculate_identity(record: dict, target: str) -> dict:
    query_position = 0
    reference_position = int(record["source_position"]) - 1
    aligned = matches = splice_junctions = 0
    for length_text, operation in CIGAR_PATTERN.findall(record["cigar"]):
        length = int(length_text)
        if operation in {"M", "=", "X"}:
            query = record["read_sequence"][query_position : query_position + length]
            reference = target[reference_position : reference_position + length]
            if len(query) != length or len(reference) != length:
                return {"coordinate_valid": False, "aligned_bases": aligned, "matches": matches, "reference_identity": None, "splice_junctions": splice_junctions, "target_end": reference_position}
            matches += sum(left == right for left, right in zip(query, reference))
            aligned += length
            query_position += length
            reference_position += length
        elif operation in {"I", "S"}:
            query_position += length
        elif operation in {"D", "N"}:
            reference_position += length
            splice_junctions += int(operation == "N")
    return {"coordinate_valid": aligned > 0, "aligned_bases": aligned, "matches": matches, "reference_identity": round(matches / aligned, 6) if aligned else None, "splice_junctions": splice_junctions, "target_end": reference_position}


def verify_protected_run(repo_root: Path, relative: str) -> int:
    run_dir = repo_root / relative
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    for name, identity in manifest["artifacts"].items():
        path = run_dir / name
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"protected evidence identity mismatch: {relative}:{name}")
    return len(manifest["artifacts"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--toolchain", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"output already exists: {args.output}")
    run = json.loads((args.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    if run["source_manifest_sha256"] != sha256(args.sources) or run["toolchain_sha256"] != sha256(args.toolchain) or run["reference_manifest_sha256"] != sha256(args.reference_manifest):
        raise SystemExit("Slice 0B input manifest identity mismatch")
    verified_artifacts = {}
    for name, identity in run["artifacts"].items():
        path = args.run_dir / name
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"Slice 0B artifact identity mismatch: {name}")
        verified_artifacts[name] = identity["sha256"]
    sources = json.loads(args.sources.read_text(encoding="utf-8"))
    acquisition_path = args.repo_root / sources["acquisition_evidence"]
    if not acquisition_path.is_file() or run.get("acquisition_evidence_sha256") != sha256(acquisition_path):
        raise SystemExit("Slice 0B acquisition evidence identity mismatch")
    source = sources["sources"][0]
    for role, identity in source["files"].items():
        path = args.repo_root / identity["path"]
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"transcript source identity mismatch: {role}")
    toolchain = json.loads(args.toolchain.read_text(encoding="utf-8"))
    for role in ("archive", "sam_dump"):
        identity = toolchain["tools"]["sra_toolkit"][role]
        path = args.repo_root / identity["path"]
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256(path) != identity["sha256"]:
            raise SystemExit(f"tool identity mismatch: {role}")
    reference_manifest = json.loads(args.reference_manifest.read_text(encoding="utf-8"))
    reference = reference_manifest["references"]["strain-b"]
    fasta_identity = reference["files"]["fasta"]
    fasta_path = args.repo_root / "reference/data/strain-b" / fasta_identity["local_name"]
    if fasta_path.stat().st_size != fasta_identity["size_bytes"] or sha256(fasta_path) != fasta_identity["sha256"]:
        raise SystemExit("primary FASTA identity mismatch")
    sequences = read_fasta(fasta_path)
    probe = json.loads((args.run_dir / "transcript_evidence_probe.json").read_text(encoding="utf-8"))
    parameters = probe["parameters"]
    if probe.get("probe_mode") != PROBE_MODE or parameters.get("probe_mode") != PROBE_MODE:
        raise SystemExit("Slice 0B probe mode mismatch")
    if parameters.get("source_sequence_map") != SEQUENCE_MAP or parameters.get("windows_per_chromosome") != list(WINDOWS):
        raise SystemExit("Slice 0B locked chromosome/window design mismatch")
    threshold = float(parameters["minimum_reference_identity"])
    endpoint_distance = int(parameters["endpoint_distance_bases"])
    recalculated_counts = Counter()
    compatible_chromosomes = set()
    endpoint_side_counts = Counter()
    endpoint_status_counts = Counter()
    endpoint_affected_records = 0
    endpoint_and_multi_mapping = 0
    for record in probe["records"]:
        target_seqid = record["target_seqid"]
        calculated = recalculate_identity(record, sequences[target_seqid])
        for key, value in calculated.items():
            if record[key] != value:
                raise SystemExit(f"record identity mismatch: {record['read_id']}:{key}")
        endpoint = reference["sequence_endpoints"][target_seqid]
        endpoint_side = endpoint_status = None
        if record["target_start"] < endpoint_distance and endpoint["five_prime"]["status"] != "complete":
            endpoint_side = "five_prime"
            endpoint_status = endpoint["five_prime"]["status"]
        elif len(sequences[target_seqid]) - record["target_end"] < endpoint_distance and endpoint["three_prime"]["status"] != "complete":
            endpoint_side = "three_prime"
            endpoint_status = endpoint["three_prime"]["status"]
        if record.get("endpoint_side") != endpoint_side or record.get("endpoint_status") != endpoint_status:
            raise SystemExit(f"record endpoint mismatch: {record['read_id']}")
        if record["nh"] > 1 or record["flag"] & 256:
            classification = "multi_mapping"
        elif not calculated["coordinate_valid"] or calculated["reference_identity"] < threshold:
            classification = "conflict"
        elif endpoint_status is not None:
            classification = "endpoint_limited"
        else:
            classification = "success"
        if classification != record["classification"]:
            raise SystemExit(f"record classification mismatch: {record['read_id']}")
        recalculated_counts[classification] += 1
        if calculated["coordinate_valid"] and calculated["reference_identity"] is not None and calculated["reference_identity"] >= threshold:
            compatible_chromosomes.add(target_seqid)
        if endpoint_status is not None:
            endpoint_affected_records += 1
            endpoint_side_counts[endpoint_side] += 1
            endpoint_status_counts[endpoint_status] += 1
            endpoint_and_multi_mapping += int(classification == "multi_mapping")
    for category in REQUIRED_CLASSES:
        recalculated_counts.setdefault(category, 0)
    summary = probe["summary"]
    if dict(sorted(recalculated_counts.items())) != summary["sample_classification_counts"] or summary["reported_category_counts"] != summary["sample_classification_counts"]:
        raise SystemExit("Slice 0B classification summary mismatch")
    expected_chromosomes = sorted(SEQUENCE_MAP.values())
    if summary["chromosomes_covered"] != expected_chromosomes or summary["coordinate_compatible_chromosomes"] != expected_chromosomes:
        raise SystemExit("Slice 0B chromosome compatibility mismatch")
    if summary.get("chromosome_coordinate_probe") != "passed" or summary.get("window_coverage") != "partial" or summary.get("coordinate_compatibility") != "partial":
        raise SystemExit("Slice 0B coordinate status mismatch")
    expected_windows = {f"{seqid}:{window}" for seqid in expected_chromosomes for window in WINDOWS}
    if set(summary["window_statistics"]) != expected_windows or summary["window_count"] != 12:
        raise SystemExit("Slice 0B window coverage mismatch")
    zero_windows = sorted(key for key, value in summary["window_statistics"].items() if value["sampled_records"] == 0)
    if zero_windows != summary["zero_coverage_windows"]:
        raise SystemExit("Slice 0B zero-coverage window mismatch")
    if summary.get("endpoint_affected_records") != endpoint_affected_records or summary.get("endpoint_and_multi_mapping") != endpoint_and_multi_mapping:
        raise SystemExit("Slice 0B endpoint count mismatch")
    if summary.get("endpoint_side_counts") != dict(sorted(endpoint_side_counts.items())) or summary.get("endpoint_status_counts") != dict(sorted(endpoint_status_counts.items())):
        raise SystemExit("Slice 0B endpoint breakdown mismatch")
    for key, value in summary["window_statistics"].items():
        records = [record for record in probe["records"] if f"{record['target_seqid']}:{record['window']}" == key]
        if value.get("endpoint_affected_records") != sum(record.get("endpoint_status") is not None for record in records) or value.get("endpoint_and_multi_mapping") != sum(record.get("endpoint_status") is not None and record["classification"] == "multi_mapping" for record in records):
            raise SystemExit(f"Slice 0B window endpoint summary mismatch: {key}")
    unavailable_fields = ("transcription_strand_status", "controlled_unique_mapping_status", "controlled_splice_evidence_status", "true_unmappable_status", "contamination_status")
    if any(not summary.get(field, "").startswith("unavailable") for field in unavailable_fields):
        raise SystemExit("Slice 0B quality limitation was overclaimed")
    diagnostics = summary.get("submitted_alignment_diagnostics", {})
    if diagnostics.get("interpretation") != "submitted NH/NM/CIGAR values are diagnostic only and are not an independent quality gate":
        raise SystemExit("Slice 0B submitted-alignment diagnostic label mismatch")
    metadata_root = ET.parse(args.repo_root / source["files"]["sra_metadata_xml"]["path"]).getroot()
    table_counts = {table.get("name", ""): int(table.find("./Statistics/Rows").get("count", "0")) for table in metadata_root.findall(".//Database/Table") if table.find("./Statistics/Rows") is not None}
    total_reads = source["library"]["read_count_per_spot"] * source["library"]["total_spots"]
    archive_stats = {"total_reads": total_reads, "primary_alignment_rows": table_counts.get("PRIMARY_ALIGNMENT", 0), "secondary_alignment_rows": table_counts.get("SECONDARY_ALIGNMENT", 0), "unmappable_or_not_primary_rows_estimate": max(0, total_reads - table_counts.get("PRIMARY_ALIGNMENT", 0))}
    if archive_stats != summary["archive_alignment_statistics"]:
        raise SystemExit("Slice 0B archive alignment statistics mismatch")
    for source_seqid, target_seqid in SEQUENCE_MAP.items():
        header = probe["reference_header_identity"][source_seqid]
        if header["normalized_length"] != len(sequences[target_seqid]) or header["raw_sam_length"] % (2**31) != len(sequences[target_seqid]):
            raise SystemExit(f"Slice 0B reference header mismatch: {source_seqid}")
    matrix = json.loads((args.run_dir / "source_applicability_matrix.json").read_text(encoding="utf-8"))["sources"][0]
    unavailable_matrix = ("controlled_unique_mapping", "controlled_splice_evidence", "true_unmappable", "contamination")
    if matrix["eligible_as_independent_transcript_source"] is not False or matrix["independence_from_current_annotation"] != "not-passed" or matrix["license_gate"] != "unresolved" or matrix["exact_target_strain_coordinates"] is not False or matrix.get("probe_mode") != PROBE_MODE or any(matrix.get(field) != "unavailable" for field in unavailable_matrix):
        raise SystemExit("Slice 0B applicability matrix overclaims source eligibility")
    if run["primary_coordinate_space"] != "GCA_001746955.1" or run["exact_target_strain_coordinates"] is not False or run["scientific_acceptance_status"] != "blocked":
        raise SystemExit("Slice 0B run scope/status mismatch")
    protected = {
        "slice0_artifact_count": verify_protected_run(args.repo_root, "local_runs/strain-b_slice0_completion_v6_run1"),
        "slice0a_artifact_count": verify_protected_run(args.repo_root, "local_runs/slice0a/qualification_v5_run1"),
    }
    result = {
        "schema_version": 1,
        "verification_type": "slice0b_independent_summary",
        "status": "passed",
        "run_id": run["run_id"],
        "verified_artifacts": verified_artifacts,
        "acquisition_evidence_sha256": sha256(acquisition_path),
        "sampled_record_count": len(probe["records"]),
        "sample_classification_counts": dict(sorted(recalculated_counts.items())),
        "chromosomes_covered": expected_chromosomes,
        "window_count": 12,
        "zero_coverage_windows": zero_windows,
        "probe_mode": PROBE_MODE,
        "chromosome_coordinate_probe": "passed",
        "window_coverage": "partial",
        "coordinate_compatibility": "partial",
        "endpoint_affected_records": endpoint_affected_records,
        "endpoint_and_multi_mapping": endpoint_and_multi_mapping,
        "quality_gate_status": "unavailable_without_controlled_read_remapping",
        "archive_alignment_statistics": archive_stats,
        "source_independence_gate": "not-passed",
        "source_license_gate": "unresolved",
        "scientific_acceptance_status": "blocked",
        "protected_evidence": protected,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
