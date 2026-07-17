from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file, write_json, write_tsv
from .mapping_probe import read_fasta_sequences
from .pipeline import _implementation_hash
from .transcript_qualification import load_transcript_source_manifest, qualify_transcript_source


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
REQUIRED_CLASSES = ("success", "conflict", "multi_mapping", "unmappable", "endpoint_limited")
EXPECTED_SEQUENCE_MAP = {
    "chr1": "CP014715.1",
    "chr2": "CP014716.1",
    "chr3": "CP014717.1",
    "chr4": "CP014718.1",
}
EXPECTED_WINDOWS = ("start", "middle", "end")
PROBE_MODE = "submitted_alignment_coordinate_probe"


def _verify_file(identity: dict[str, Any], repo_root: Path, label: str) -> Path:
    path = repo_root / identity["path"]
    if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256_file(path) != identity["sha256"]:
        raise ContractError(f"Slice 0B file identity mismatch: {label}")
    return path


def _sam_dump_records(executable: Path, archive: Path, region: str, limit: int) -> list[str]:
    process = subprocess.Popen(
        [str(executable), "--no-header", "--aligned-region", region, str(archive)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None
    records: list[str] = []
    try:
        for raw in process.stdout:
            line = raw.rstrip("\r\n")
            if line and not line.startswith("@"):
                records.append(line)
                if len(records) >= limit:
                    process.terminate()
                    break
        process.wait(timeout=30)
    except Exception:
        process.kill()
        process.wait(timeout=10)
        raise
    if len(records) < limit and process.returncode not in (0, None):
        raise ContractError(f"sam-dump failed for region {region}: exit {process.returncode}")
    return records


def _sam_header(executable: Path, archive: Path) -> dict[str, int]:
    completed = subprocess.run(
        [str(executable), "--aligned-region", "chr1:1-1", str(archive)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    if completed.returncode != 0:
        raise ContractError("sam-dump header probe failed")
    result = {}
    for line in completed.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(item.split(":", 1) for item in line.split("\t")[1:])
        result[fields["SN"]] = int(fields["LN"])
    if not result:
        raise ContractError("sam-dump header has no reference identities")
    return result


def _parse_sam(line: str) -> dict[str, Any]:
    columns = line.split("\t")
    if len(columns) < 11:
        raise ContractError("malformed SAM record")
    tags = {}
    for item in columns[11:]:
        parts = item.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return {
        "read_id": columns[0],
        "flag": int(columns[1]),
        "source_seqid": columns[2],
        "source_position": int(columns[3]),
        "mapq": int(columns[4]),
        "cigar": columns[5],
        "mate_seqid": columns[6],
        "mate_position": int(columns[7]),
        "template_length": int(columns[8]),
        "read_sequence": columns[9],
        "nh": int(tags.get("NH", "1")),
        "nm": int(tags["NM"]) if "NM" in tags else None,
    }


def _reference_identity(record: dict[str, Any], target_sequence: str) -> dict[str, Any]:
    query_position = 0
    reference_position = record["source_position"] - 1
    aligned = matches = splice_junctions = 0
    for length_text, operation in CIGAR_PATTERN.findall(record["cigar"]):
        length = int(length_text)
        if operation in {"M", "=", "X"}:
            query = record["read_sequence"][query_position : query_position + length]
            reference = target_sequence[reference_position : reference_position + length]
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
            splice_junctions += operation == "N"
        elif operation in {"H", "P"}:
            continue
    return {
        "coordinate_valid": aligned > 0,
        "aligned_bases": aligned,
        "matches": matches,
        "reference_identity": round(matches / aligned, 6) if aligned else None,
        "splice_junctions": splice_junctions,
        "target_end": reference_position,
    }


def _archive_alignment_statistics(metadata_path: Path, read_count_per_spot: int, total_spots: int) -> dict[str, int]:
    root = ET.parse(metadata_path).getroot()
    table_counts = {
        table.get("name", ""): int(table.find("./Statistics/Rows").get("count", "0"))
        for table in root.findall(".//Database/Table")
        if table.find("./Statistics/Rows") is not None
    }
    total_reads = read_count_per_spot * total_spots
    primary = table_counts.get("PRIMARY_ALIGNMENT", 0)
    return {
        "total_reads": total_reads,
        "primary_alignment_rows": primary,
        "secondary_alignment_rows": table_counts.get("SECONDARY_ALIGNMENT", 0),
        "unmappable_or_not_primary_rows_estimate": max(0, total_reads - primary),
    }


def run_transcript_probe(
    source_manifest_path: Path,
    toolchain_path: Path,
    reference_manifest_path: Path,
    output_dir: Path,
    repo_root: Path,
) -> dict[str, Any]:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise ContractError(f"output directory already exists: {output_dir}")
    source_manifest = load_transcript_source_manifest(source_manifest_path)
    if len(source_manifest["sources"]) != 1:
        raise ContractError("current Slice 0B probe requires exactly one qualified source")
    source = source_manifest["sources"][0]
    source_qualification = qualify_transcript_source(source, repo_root)
    toolchain = json.loads(toolchain_path.read_text(encoding="utf-8"))
    tool = toolchain["tools"]["sra_toolkit"]
    _verify_file(tool["archive"], repo_root, "SRA Toolkit archive")
    sam_dump = _verify_file(tool["sam_dump"], repo_root, "sam-dump")
    version = subprocess.run([str(sam_dump), "--version"], capture_output=True, text=True, encoding="utf-8", check=False)
    if version.returncode != 0 or tool["version"] not in version.stdout:
        raise ContractError("sam-dump version mismatch")
    archive = repo_root / source["files"]["sra_lite"]["path"]
    reference_manifest = json.loads(reference_manifest_path.read_text(encoding="utf-8"))
    reference = reference_manifest["references"]["strain-b"]
    if reference["assembly_accession"] != "GCA_001746955.1":
        raise ContractError("Slice 0B probe changed the locked primary assembly")
    fasta_identity = reference["files"]["fasta"]
    fasta_path = repo_root / "reference/data/strain-b" / fasta_identity["local_name"]
    if not fasta_path.is_file() or fasta_path.stat().st_size != fasta_identity["size_bytes"] or sha256_file(fasta_path) != fasta_identity["sha256"]:
        raise ContractError("Slice 0B primary FASTA identity mismatch")
    sequences = read_fasta_sequences(fasta_path)
    parameters = toolchain["probe_parameters"]
    sequence_map = parameters["source_sequence_map"]
    if parameters.get("probe_mode") != PROBE_MODE:
        raise ContractError("Slice 0B probe mode must remain a submitted-alignment coordinate probe")
    if sequence_map != EXPECTED_SEQUENCE_MAP or parameters.get("windows_per_chromosome") != list(EXPECTED_WINDOWS):
        raise ContractError("Slice 0B probe must use the locked four-chromosome and 12-window design")
    header = _sam_header(sam_dump, archive)
    header_checks = {}
    for source_seqid, target_seqid in sequence_map.items():
        normalized_length = header[source_seqid] % (2**31)
        target_length = len(sequences[target_seqid])
        header_checks[source_seqid] = {
            "raw_sam_length": header[source_seqid],
            "normalized_length": normalized_length,
            "target_seqid": target_seqid,
            "target_length": target_length,
            "length_identity": normalized_length == target_length,
        }
    if not all(item["length_identity"] for item in header_checks.values()):
        raise ContractError("SRA alignment reference lengths do not match GCA_001746955.1")
    window_length = int(parameters["window_length_bases"])
    limit = int(parameters["max_sam_records_per_window"])
    identity_threshold = float(parameters["minimum_reference_identity"])
    endpoint_distance = int(parameters["endpoint_distance_bases"])
    records: list[dict[str, Any]] = []
    designed_windows: dict[tuple[str, str], dict[str, Any]] = {}
    for source_seqid, target_seqid in sequence_map.items():
        target = sequences[target_seqid]
        length = len(target)
        windows = {
            "start": (1, min(length, window_length)),
            "middle": (max(1, length // 2 - window_length // 2 + 1), min(length, length // 2 + window_length // 2)),
            "end": (max(1, length - window_length + 1), length),
        }
        endpoint_record = reference["sequence_endpoints"][target_seqid]
        for window_name, (window_start, window_end) in windows.items():
            region = f"{source_seqid}:{window_start}-{window_end}"
            designed_windows[(target_seqid, window_name)] = {
                "source_seqid": source_seqid,
                "region": region,
                "window_start": window_start,
                "window_end": window_end,
            }
            for line in _sam_dump_records(sam_dump, archive, region, limit):
                record = _parse_sam(line)
                identity = _reference_identity(record, target)
                record.update(identity)
                record.update({"window": window_name, "region": region, "target_seqid": target_seqid, "target_start": record["source_position"] - 1, "alignment_orientation": "reverse" if record["flag"] & 16 else "forward"})
                endpoint_side = endpoint_status = None
                if record["target_start"] < endpoint_distance and endpoint_record["five_prime"]["status"] != "complete":
                    endpoint_side = "five_prime"
                    endpoint_status = endpoint_record["five_prime"]["status"]
                elif length - record["target_end"] < endpoint_distance and endpoint_record["three_prime"]["status"] != "complete":
                    endpoint_side = "three_prime"
                    endpoint_status = endpoint_record["three_prime"]["status"]
                record["endpoint_side"] = endpoint_side
                record["endpoint_status"] = endpoint_status
                if record["nh"] > 1 or record["flag"] & 256:
                    category = "multi_mapping"
                elif not record["coordinate_valid"] or record["reference_identity"] < identity_threshold:
                    category = "conflict"
                elif endpoint_status is not None:
                    category = "endpoint_limited"
                else:
                    category = "success"
                record["classification"] = category
                records.append(record)
    records.sort(key=lambda item: (item["target_seqid"], item["window"], item["source_position"], item["read_id"], item["flag"]))
    sample_counts = Counter(item["classification"] for item in records)
    for category in REQUIRED_CLASSES:
        sample_counts.setdefault(category, 0)
    metadata_path = repo_root / source["files"]["sra_metadata_xml"]["path"]
    archive_statistics = _archive_alignment_statistics(metadata_path, source["library"]["read_count_per_spot"], source["library"]["total_spots"])
    reported_counts = dict(sorted(sample_counts.items()))
    by_window: dict[str, dict[str, Any]] = {}
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[(record["target_seqid"], record["window"])].append(record)
    for target_seqid, window in sorted(designed_windows):
        items = grouped[(target_seqid, window)]
        identities = [item["reference_identity"] for item in items if item["reference_identity"] is not None]
        by_window[f"{target_seqid}:{window}"] = {
            **designed_windows[(target_seqid, window)],
            "sampled_records": len(items),
            "classification_counts": dict(sorted(Counter(item["classification"] for item in items).items())),
            "median_reference_identity": round(median(identities), 6) if identities else None,
            "splice_junction_records": sum(item["splice_junctions"] > 0 for item in items),
            "endpoint_affected_records": sum(item["endpoint_status"] is not None for item in items),
            "endpoint_and_multi_mapping": sum(item["endpoint_status"] is not None and item["classification"] == "multi_mapping" for item in items),
        }
    chromosome_coverage = sorted({item["target_seqid"] for item in records})
    compatible_chromosomes = sorted({
        item["target_seqid"]
        for item in records
        if item["coordinate_valid"] and item["reference_identity"] is not None and item["reference_identity"] >= identity_threshold
    })
    chromosome_coordinate_probe = "passed" if compatible_chromosomes == sorted(EXPECTED_SEQUENCE_MAP.values()) else "failed"
    zero_coverage_windows = sorted(key for key, value in by_window.items() if value["sampled_records"] == 0)
    windows_with_records = len(by_window) - len(zero_coverage_windows)
    window_coverage = "passed" if windows_with_records == 12 else "partial" if windows_with_records > 0 else "failed"
    coordinate_compatibility = "passed" if chromosome_coordinate_probe == "passed" and window_coverage == "passed" else "partial" if chromosome_coordinate_probe == "passed" and window_coverage == "partial" else "failed"
    endpoint_records = [item for item in records if item["endpoint_status"] is not None]
    probe = {
        "schema_version": 1,
        "source_id": source["source_id"],
        "source_archive": source["files"]["sra_lite"],
        "primary_coordinate_space": "GCA_001746955.1",
        "exact_target_strain_coordinates": False,
        "probe_mode": PROBE_MODE,
        "tool": {"name": tool["name"], "version": tool["version"], "sam_dump": tool["sam_dump"]},
        "parameters": parameters,
        "reference_header_identity": header_checks,
        "records": records,
        "summary": {
            "sampled_record_count": len(records),
            "sample_classification_counts": dict(sorted(sample_counts.items())),
            "reported_category_counts": reported_counts,
            "archive_alignment_statistics": archive_statistics,
            "true_unmappable_status": "unavailable; submitted alignments cannot independently establish true unmappable reads",
            "chromosomes_covered": chromosome_coverage,
            "coordinate_compatible_chromosomes": compatible_chromosomes,
            "window_count": len(by_window),
            "windows_with_records": windows_with_records,
            "zero_coverage_windows": zero_coverage_windows,
            "window_statistics": by_window,
            "endpoint_affected_records": len(endpoint_records),
            "endpoint_and_multi_mapping": sum(item["classification"] == "multi_mapping" for item in endpoint_records),
            "endpoint_side_counts": dict(sorted(Counter(item["endpoint_side"] for item in endpoint_records).items())),
            "endpoint_status_counts": dict(sorted(Counter(item["endpoint_status"] for item in endpoint_records).items())),
            "submitted_alignment_diagnostics": {
                "nh_unique_records": sum(item["nh"] == 1 and not item["flag"] & 256 for item in records),
                "multi_mapping_records": sum(item["classification"] == "multi_mapping" for item in records),
                "splice_junction_records": sum(item["splice_junctions"] > 0 for item in records),
                "alignment_orientation_counts": dict(sorted(Counter(item["alignment_orientation"] for item in records).items())),
                "interpretation": "submitted NH/NM/CIGAR values are diagnostic only and are not an independent quality gate",
            },
            "transcription_strand_status": "unavailable; library is stranded but first/second-strand orientation convention is not stated",
            "controlled_unique_mapping_status": "unavailable; no controlled read remapping was performed",
            "controlled_splice_evidence_status": "unavailable; submitted CIGAR N operations are diagnostic only",
            "contamination_status": "unavailable; submitted alignments cannot independently establish contamination",
            "chromosome_coordinate_probe": chromosome_coordinate_probe,
            "window_coverage": window_coverage,
            "coordinate_compatibility": coordinate_compatibility,
        },
    }
    applicability = [{
        "source_id": source["source_id"],
        "archive_strain": source["strain"],
        "publication_and_assembly_strain": "Strain-B, ATCC 20864",
        "target_strain": "Strain-T",
        "exact_target_strain_coordinates": False,
        "strain_secondary_identifier_status": "conflict",
        "annotation_provenance_overlap": True,
        "independence_from_current_annotation": "not-passed",
        "file_license_assignment": source["license"]["file_license_assignment"],
        "license_gate": "unresolved",
        "strand_specificity": "passed",
        "strand_orientation": "unavailable",
        "probe_mode": PROBE_MODE,
        "controlled_unique_mapping": "unavailable",
        "controlled_splice_evidence": "unavailable",
        "true_unmappable": "unavailable",
        "contamination": "unavailable",
        "coordinate_compatibility": probe["summary"]["coordinate_compatibility"],
        "eligible_as_independent_transcript_source": False,
    }]
    implementation_hash = _implementation_hash()
    material = json.dumps({
        "source_manifest_sha256": sha256_file(source_manifest_path),
        "toolchain_sha256": sha256_file(toolchain_path),
        "reference_manifest_sha256": sha256_file(reference_manifest_path),
        "implementation_sha256": implementation_hash,
        "parameters": parameters,
    }, sort_keys=True).encode("utf-8")
    run_id = "slice0b-" + hashlib.sha256(material).hexdigest()[:16]
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        write_json(temp_dir / "transcript_evidence_probe.json", probe)
        fields = sorted({key for record in records for key in record})
        write_tsv(temp_dir / "transcript_evidence_probe.tsv", records, fields)
        write_json(temp_dir / "source_applicability_matrix.json", {"schema_version": 1, "sources": applicability})
        write_tsv(temp_dir / "source_applicability_matrix.tsv", applicability, list(applicability[0]))
        report = "\n".join([
            "# Slice 0B transcript evidence probe",
            "",
            "## Scope",
            "",
            "- Primary coordinate space remains Strain-B GCA_001746955.1.",
            "- Coordinates are not exact Strain-T coordinates.",
            "- No safe-harbor candidates or thresholds are generated.",
            "",
            "## Probe result",
            "",
            f"- Sampled SAM records: {len(records)} across {len(chromosome_coverage)} nuclear chromosomes; {len(by_window)} windows were designed, {len(by_window) - len(zero_coverage_windows)} contained records and {len(zero_coverage_windows)} had zero sampled coverage.",
            f"- Sampled categories: {json.dumps(reported_counts, sort_keys=True)}.",
            f"- Archive unmappable-or-not-primary estimate: {archive_statistics['unmappable_or_not_primary_rows_estimate']}; this is not a verified per-read unmappable count.",
            f"- Endpoint-affected records: {probe['summary']['endpoint_affected_records']}; endpoint plus multi-mapping records: {probe['summary']['endpoint_and_multi_mapping']}. The mutually exclusive endpoint_limited category alone does not represent all endpoint risk.",
            f"- Submitted-alignment diagnostics: {json.dumps(probe['summary']['submitted_alignment_diagnostics'], sort_keys=True)}.",
            f"- Chromosome coordinate probe: {probe['summary']['chromosome_coordinate_probe']}; window coverage: {probe['summary']['window_coverage']}; overall coordinate compatibility: {probe['summary']['coordinate_compatibility']}.",
            f"- Zero-coverage windows: {json.dumps(zero_coverage_windows)}.",
            "- Transcription strand remains unavailable because the archive does not state the orientation convention.",
            "- Unmapped/non-primary rows are not interpreted as contamination.",
            "- No controlled read remapping was performed; controlled uniqueness, splice evidence, true unmappable reads and contamination remain unavailable.",
            "",
            "## Qualification conclusion",
            "",
            "The run is a useful coordinate-positive stranded RNA-seq control, but it is not an independent source for the current annotation: its transcript models contributed to initial assembly/annotation. The ATCC identifier conflicts and the per-file license remains unresolved. Scientific status therefore remains blocked.",
        ]) + "\n"
        (temp_dir / "transcript_quality_report.md").write_text(report, encoding="utf-8")
        recommendation = "\n".join([
            "# Slice 0B recommendation",
            "",
            "Recommendation: evidence is insufficient to authorize a full-genome boundary track.",
            "",
            "Reasons:",
            "",
            "- PRJNA311606 provides a strong same-study coordinate control but overlaps the provenance of the current assembly and initial annotation.",
            "- The archive and publication disagree on the ATCC secondary strain identifier.",
            "- No explicit per-file SRA license was found; NCBI database usage terms are not promoted to a file license.",
            "- The library is stranded, but the orientation convention is unavailable.",
            "- The coordinate probe uses submitted alignments; controlled uniqueness, splice evidence, true unmappable reads and contamination remain unavailable.",
            "",
            "Scientific status: blocked. Do not start full-genome reannotation or Slice 1.",
        ]) + "\n"
        (temp_dir / "slice0b_recommendation.md").write_text(recommendation, encoding="utf-8")
        artifacts = [
            "transcript_evidence_probe.json",
            "transcript_evidence_probe.tsv",
            "source_applicability_matrix.json",
            "source_applicability_matrix.tsv",
            "transcript_quality_report.md",
            "slice0b_recommendation.md",
        ]
        run_manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "slice": "Slice 0B experimental evidence qualification",
            "execution_status": "complete",
            "verification_status": "not_run",
            "scientific_acceptance_status": "blocked",
            "scientific_acceptance_blockers": [
                "no transcript source independent of the current annotation-generation process has passed",
                "per-file transcript data license unresolved",
                "archive/publication ATCC identifier conflict",
            ],
            "primary_coordinate_space": "GCA_001746955.1",
            "exact_target_strain_coordinates": False,
            "source_manifest_sha256": sha256_file(source_manifest_path),
            "acquisition_evidence_sha256": sha256_file(repo_root / source_manifest["acquisition_evidence"]),
            "toolchain_sha256": sha256_file(toolchain_path),
            "reference_manifest_sha256": sha256_file(reference_manifest_path),
            "implementation_sha256": implementation_hash,
            "artifacts": {name: {"size_bytes": (temp_dir / name).stat().st_size, "sha256": sha256_file(temp_dir / name)} for name in artifacts},
        }
        write_json(temp_dir / "run_manifest.json", run_manifest)
        os.replace(temp_dir, output_dir)
        return run_manifest
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


def create_slice0b_acceptance(
    run_dir: Path,
    independent_path: Path,
    test_evidence_path: Path,
    repeatability_path: Path,
    peer_run_dir: Path,
    repo_root: Path,
    source_manifest_path: Path,
    toolchain_path: Path,
    reference_manifest_path: Path,
) -> dict[str, Any]:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    independent = json.loads(independent_path.read_text(encoding="utf-8"))
    tests = json.loads(test_evidence_path.read_text(encoding="utf-8"))
    repeatability = json.loads(repeatability_path.read_text(encoding="utf-8"))
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    expected_artifacts = {
        "transcript_evidence_probe.json",
        "transcript_evidence_probe.tsv",
        "source_applicability_matrix.json",
        "source_applicability_matrix.tsv",
        "transcript_quality_report.md",
        "slice0b_recommendation.md",
    }
    if set(run.get("artifacts", {})) != expected_artifacts:
        raise ContractError("Slice 0B artifact set mismatch")
    if run.get("primary_coordinate_space") != "GCA_001746955.1" or run.get("exact_target_strain_coordinates") is not False:
        raise ContractError("Slice 0B coordinate scope mismatch")
    if run.get("execution_status") != "complete" or run.get("scientific_acceptance_status") != "blocked":
        raise ContractError("Slice 0B run status mismatch")
    if run.get("implementation_sha256") != _implementation_hash():
        raise ContractError("Slice 0B implementation identity is stale")
    if run.get("source_manifest_sha256") != sha256_file(source_manifest_path) or run.get("toolchain_sha256") != sha256_file(toolchain_path) or run.get("reference_manifest_sha256") != sha256_file(reference_manifest_path):
        raise ContractError("Slice 0B input manifest identity mismatch")
    acquisition_relative = source_manifest.get("acquisition_evidence")
    if not isinstance(acquisition_relative, str) or not acquisition_relative:
        raise ContractError("Slice 0B acquisition evidence path is missing")
    resolved_repo_root = repo_root.resolve()
    acquisition_path = (resolved_repo_root / acquisition_relative).resolve()
    if not acquisition_path.is_relative_to(resolved_repo_root):
        raise ContractError("Slice 0B acquisition evidence escaped the repository root")
    if not acquisition_path.is_file() or run.get("acquisition_evidence_sha256") != sha256_file(acquisition_path):
        raise ContractError("Slice 0B acquisition evidence identity mismatch")
    acquisition_identity = {"path": acquisition_relative, "size_bytes": acquisition_path.stat().st_size, "sha256": sha256_file(acquisition_path)}
    for name, identity in run["artifacts"].items():
        path = run_dir / name
        if not path.is_file() or path.stat().st_size != identity["size_bytes"] or sha256_file(path) != identity["sha256"]:
            raise ContractError(f"Slice 0B artifact identity mismatch: {name}")
    if independent.get("schema_version") != 1 or independent.get("verification_type") != "slice0b_independent_summary" or independent.get("status") != "passed" or independent.get("run_id") != run["run_id"]:
        raise ContractError("Slice 0B independent evidence identity mismatch")
    if independent.get("verified_artifacts") != {name: identity["sha256"] for name, identity in run["artifacts"].items()}:
        raise ContractError("Slice 0B independent artifact summary mismatch")
    if tests.get("schema_version") != 1 or tests.get("evidence_type") != "automated_tests" or tests.get("status") != "passed" or tests.get("command") != ["python", "-m", "pytest", "-q"] or tests.get("passed_count", 0) < 52 or tests.get("implementation_sha256") != run["implementation_sha256"]:
        raise ContractError("Slice 0B automated test evidence mismatch")
    probe = json.loads((run_dir / "transcript_evidence_probe.json").read_text(encoding="utf-8"))
    summary = probe["summary"]
    expected_chromosomes = ["CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"]
    if set(summary["reported_category_counts"]) != set(REQUIRED_CLASSES):
        raise ContractError("Slice 0B required evidence categories are missing")
    if probe.get("probe_mode") != PROBE_MODE or summary["chromosomes_covered"] != expected_chromosomes or summary["coordinate_compatible_chromosomes"] != expected_chromosomes:
        raise ContractError("Slice 0B four-chromosome coordinate gate mismatch")
    if summary.get("chromosome_coordinate_probe") != "passed" or summary.get("window_coverage") != "partial" or summary.get("coordinate_compatibility") != "partial":
        raise ContractError("Slice 0B coordinate status was overclaimed")
    if summary["window_count"] != 12 or set(summary["window_statistics"]) != {f"{seqid}:{window}" for seqid in expected_chromosomes for window in ("start", "middle", "end")}:
        raise ContractError("Slice 0B window coverage mismatch")
    unavailable_fields = ("transcription_strand_status", "controlled_unique_mapping_status", "controlled_splice_evidence_status", "true_unmappable_status", "contamination_status")
    if any(not summary.get(field, "").startswith("unavailable") for field in unavailable_fields):
        raise ContractError("Slice 0B quality limitation report mismatch")
    endpoint_records = [item for item in probe["records"] if item.get("endpoint_status") is not None]
    if summary.get("endpoint_affected_records") != len(endpoint_records) or summary.get("endpoint_and_multi_mapping") != sum(item["classification"] == "multi_mapping" for item in endpoint_records):
        raise ContractError("Slice 0B endpoint risk summary mismatch")
    matrix = json.loads((run_dir / "source_applicability_matrix.json").read_text(encoding="utf-8"))["sources"]
    if len(matrix) != 1:
        raise ContractError("Slice 0B applicability matrix source count mismatch")
    source = matrix[0]
    unavailable_matrix = ("controlled_unique_mapping", "controlled_splice_evidence", "true_unmappable", "contamination")
    if source.get("eligible_as_independent_transcript_source") is not False or source.get("independence_from_current_annotation") != "not-passed" or source.get("license_gate") != "unresolved" or source.get("strain_secondary_identifier_status") != "conflict" or source.get("probe_mode") != PROBE_MODE or any(source.get(field) != "unavailable" for field in unavailable_matrix):
        raise ContractError("Slice 0B source eligibility was overclaimed")
    if independent.get("source_independence_gate") != "not-passed" or independent.get("source_license_gate") != "unresolved" or independent.get("scientific_acceptance_status") != "blocked":
        raise ContractError("Slice 0B independent scientific gate mismatch")
    independent_expected = {
        "probe_mode": PROBE_MODE,
        "sampled_record_count": summary["sampled_record_count"],
        "sample_classification_counts": summary["sample_classification_counts"],
        "chromosomes_covered": expected_chromosomes,
        "window_count": 12,
        "zero_coverage_windows": summary["zero_coverage_windows"],
        "chromosome_coordinate_probe": "passed",
        "window_coverage": "partial",
        "coordinate_compatibility": "partial",
        "endpoint_affected_records": summary["endpoint_affected_records"],
        "endpoint_and_multi_mapping": summary["endpoint_and_multi_mapping"],
        "quality_gate_status": "unavailable_without_controlled_read_remapping",
        "archive_alignment_statistics": summary["archive_alignment_statistics"],
        "acquisition_evidence_sha256": acquisition_identity["sha256"],
    }
    if any(independent.get(key) != value for key, value in independent_expected.items()):
        raise ContractError("Slice 0B independent recomputation summary mismatch")
    if independent.get("protected_evidence") != {"slice0_artifact_count": 9, "slice0a_artifact_count": 7}:
        raise ContractError("Slice 0B protected evidence check mismatch")
    repeat_files = {
        *expected_artifacts,
        "run_manifest.json",
        "verification/independent_check.json",
        "verification/test_evidence.json",
    }
    if repeatability.get("schema_version") != 1 or repeatability.get("verification_type") != "slice0b_repeatability" or repeatability.get("status") != "passed" or repeatability.get("run_id") != run["run_id"] or repeatability.get("all_files_identical") is not True:
        raise ContractError("Slice 0B repeatability evidence mismatch")
    run_directories = set(repeatability.get("run_directories", []))
    if run_dir.resolve() == peer_run_dir.resolve() or run_directories != {run_dir.name, peer_run_dir.name}:
        raise ContractError("Slice 0B repeatability run identity mismatch")
    if set(repeatability.get("file_sha256", {})) != repeat_files:
        raise ContractError("Slice 0B repeatability file set mismatch")
    peer_files = {path.relative_to(peer_run_dir).as_posix() for path in peer_run_dir.rglob("*") if path.is_file()}
    if peer_files != repeat_files:
        raise ContractError("Slice 0B peer run file set mismatch")
    for relative, digest in repeatability["file_sha256"].items():
        if sha256_file(run_dir / Path(relative)) != digest or sha256_file(peer_run_dir / Path(relative)) != digest:
            raise ContractError(f"Slice 0B repeatability hash mismatch: {relative}")
    report = (run_dir / "transcript_quality_report.md").read_text(encoding="utf-8")
    recommendation = (run_dir / "slice0b_recommendation.md").read_text(encoding="utf-8")
    if "No safe-harbor candidates or thresholds are generated" not in report or "Scientific status: blocked" not in recommendation or "Do not start full-genome reannotation or Slice 1" not in recommendation:
        raise ContractError("Slice 0B report stop line mismatch")
    result = {
        "schema_version": 1,
        "acceptance_version": "slice0b-adr0005-v1",
        "run_id": run["run_id"],
        "execution_status": "complete",
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
        "scientific_acceptance_blockers": run["scientific_acceptance_blockers"],
        "primary_coordinate_space": "GCA_001746955.1",
        "exact_target_strain_coordinates": False,
        "independent_transcript_source_passed": False,
        "acquisition_evidence": acquisition_identity,
        "run_manifest": {"size_bytes": (run_dir / "run_manifest.json").stat().st_size, "sha256": sha256_file(run_dir / "run_manifest.json")},
        "artifacts": run["artifacts"],
        "verification_evidence": {
            "independent": {"size_bytes": independent_path.stat().st_size, "sha256": sha256_file(independent_path)},
            "automated_tests": {"size_bytes": test_evidence_path.stat().st_size, "sha256": sha256_file(test_evidence_path), "passed_count": tests["passed_count"]},
            "repeatability": {"size_bytes": repeatability_path.stat().st_size, "sha256": sha256_file(repeatability_path), "compared_file_count": len(repeat_files)},
        },
        "protected_evidence": independent["protected_evidence"],
        "next_authoritative_action": "remain blocked and obtain an independently licensed transcript source that did not contribute to the current annotation-generation process",
        "stop_line": "do not start full-genome reannotation, safe-harbor candidate generation, threshold freezing, Streamlit, or Slice 1",
    }
    write_json(run_dir / "acceptance_manifest.json", result)
    return result
