from __future__ import annotations

import json
import hashlib
import subprocess
import sys
from pathlib import Path

import jsonschema

from scripts.build_srr31989028_probe_manifest import digest
from scripts.prepare_srr31989028_reads import fastq_pair_stats, forbidden_replicate_files, metadata_sra_identity
from scripts.run_srr31989028_controlled_probe import parse_summary, verified_index_prefix, wsl_path


ROOT = Path(__file__).resolve().parents[1]


def test_alignment_summary_keeps_unique_multi_unmapped_and_splice_separate(tmp_path: Path) -> None:
    sequence = list("A" * 40); sequence[5:7] = "GT"; sequence[13:15] = "AG"
    fasta = tmp_path / "ref.fa"; fasta.write_text(">CP014715.1\n" + "".join(sequence) + "\n", encoding="utf-8")
    sam = "\n".join([
        "q1\t99\tCP014715.1\t1\t60\t5M10N5M\t=\t21\t30\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:1",
        "q1\t147\tCP014715.1\t21\t60\t10M\t=\t1\t-30\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:1",
        "q2\t99\tCP014715.1\t2\t20\t10M\t=\t22\t30\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:2",
        "q2\t147\tCP014715.1\t22\t20\t10M\t=\t2\t-30\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:2",
        "q2\t355\tCP014716.1\t3\t10\t10M\t=\t23\t30\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:2",
        "q3\t77\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
        "q3\t141\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
    ]) + "\n"
    summary = tmp_path / "summary.json"; junctions = tmp_path / "junctions.tsv"
    completed = subprocess.run([sys.executable, str(ROOT / "scripts/summarize_srr31989028_alignments.py"), "--reference-fasta", str(fasta), "--summary-json", str(summary), "--junctions-tsv", str(junctions)], input=sam, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    value = json.loads(summary.read_text(encoding="utf-8")); counts = value["counts"]
    assert counts["total_query_pairs"] == 3
    assert counts["pair_unique"] == 1
    assert counts["pair_multi"] == 1
    assert counts["pair_unmapped"] == 1
    assert value["classification_complete"] is True
    assert "GT-AG\ttrue\t1\t0\tunavailable" in junctions.read_text(encoding="utf-8")


def test_coverage_summary_requires_nonzero_four_chromosome_evidence(tmp_path: Path) -> None:
    refs = ["CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"]
    config = {"sequence_lengths": {name: 10 for name in refs}, "windows": [{"window_id": f"{name}:all", "reference": name, "start_0based": 0, "end_0based": 10, "position_class": "test", "endpoint_caveat": "test"} for name in refs]}
    config_path = tmp_path / "windows.json"; config_path.write_text(json.dumps(config), encoding="utf-8")
    depth = "".join(f"{name}\t{pos}\t{1 if pos == 1 else 0}\n" for name in refs for pos in range(1, 11))
    summary = tmp_path / "coverage.json"
    completed = subprocess.run([sys.executable, str(ROOT / "scripts/summarize_srr31989028_coverage.py"), "--windows-json", str(config_path), "--summary-json", str(summary), "--windows-tsv", str(tmp_path / "windows.tsv"), "--anomalies-tsv", str(tmp_path / "anomalies.tsv")], input=depth, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(summary.read_text(encoding="utf-8"))["nuclear_chromosomes_nonzero"] is True


def test_hisat2_summary_parser_does_not_infer_strand(tmp_path: Path) -> None:
    path = tmp_path / "summary.txt"; path.write_text("Total pairs: 10\nAligned concordantly 1 time: 2\nOverall alignment rate: 20.00%\n", encoding="utf-8")
    value = parse_summary(path)
    assert value["total_pairs"] == 10
    assert value["overall_alignment_rate"] == 20.0
    assert all("strand" not in key for key in value)


def test_probe_schemas_are_valid() -> None:
    for name in ("srr31989028_acquisition_evidence.v1.schema.json", "srr31989028_probe_manifest.v1.schema.json", "srr31989028_probe_acceptance.v1.schema.json", "srr31989028_probe_failure_acceptance.v1.schema.json", "srr31989028_probe_failure_repeatability.v1.schema.json"):
        schema = json.loads((ROOT / "schemas" / name).read_text(encoding="utf-8"))
        jsonschema.Draft202012Validator.check_schema(schema)


def test_acquisition_schema_accepts_full_normalized_sra_identity() -> None:
    schema = json.loads((ROOT / "schemas/srr31989028_acquisition_evidence.v1.schema.json").read_text(encoding="utf-8"))
    value = {
        "path": "local_runs/controlled_probe/prjna1210090_srr31989028/raw/SRR31989028.sra",
        "size_bytes": 3931837205,
        "sha256": "0" * 64,
        "md5": "8c8a58880254746890e10c68b402c875",
        "format": "SRA Normalized",
        "lite": False,
    }
    jsonschema.validate(value, schema["properties"]["sra"])


def test_acquisition_schema_accepts_complete_evidence_shape() -> None:
    schema = json.loads((ROOT / "schemas/srr31989028_acquisition_evidence.v1.schema.json").read_text(encoding="utf-8"))
    identity = {"path": "local_runs/evidence.bin", "size_bytes": 1, "sha256": "0" * 64}
    command = {"argv": ["tool", "arg"], "cwd": str(ROOT), "started_at_utc": "2026-07-15T00:00:00+00:00", "completed_at_utc": "2026-07-15T00:00:01+00:00", "exit_code": 0, "stdout": identity, "stderr": identity, "status": "passed"}
    value = {
        "schema_version": 1, "acquisition_method": "ncbi-aws-normalized-sra-to-deterministic-paired-fastq", "run_accession": "SRR31989028",
        "source": {"url": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR31989028/SRR31989028", "declared_size_bytes": 3931837205, "declared_md5": "8c8a58880254746890e10c68b402c875", "metadata_snapshot": identity},
        "download_receipt": identity,
        "sra": {"path": "local_runs/controlled_probe/prjna1210090_srr31989028/raw/SRR31989028.sra", "size_bytes": 3931837205, "sha256": "1" * 64, "md5": "8c8a58880254746890e10c68b402c875", "format": "SRA Normalized", "lite": False},
        "validation": command | {"tool": identity, "version": "3.4.1"},
        "conversion": command | {"tool": identity, "version": "3.4.1", "threads": 8, "uncompressed_reads": {"read1": {"size_bytes": 1, "sha256": "2" * 64}, "read2": {"size_bytes": 1, "sha256": "3" * 64}}, "statistics": {"pair_count": 21957417, "read_count_per_file": 21957417, "total_bases": 6587225100}},
        "compression": command | {"tool": "python-standard-library:gzip.GzipFile", "version": "Python 3.12.10; zlib 1.3.1", "parameters": ["compresslevel=6", "mtime=0", "filename="]},
        "reads": {"read1": identity, "read2": identity}, "forbidden_replicates_present": False, "completed_at_utc": "2026-07-15T00:00:02+00:00",
    }
    jsonschema.validate(value, schema)


def test_wsl_path_is_bound_to_windows_drive() -> None:
    assert wsl_path(ROOT).startswith("/mnt/c/")


def test_probe_digest_streams_actual_file_content(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq.gz"; path.write_bytes(b"actual-fastq-content")
    assert digest(path, "md5") == hashlib.md5(b"actual-fastq-content").hexdigest()
    assert digest(path, "sha256") == hashlib.sha256(b"actual-fastq-content").hexdigest()


def test_ncbi_normalized_sra_metadata_identity_is_bound() -> None:
    snapshot = ROOT / "local_runs/independent_transcript_sources/prjna1210090/source_files/sra_experiment_packages.xml"
    assert metadata_sra_identity(snapshot) == {
        "url": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR31989028/SRR31989028",
        "size_bytes": 3931837205,
        "md5": "8c8a58880254746890e10c68b402c875",
    }


def test_fastq_pair_stats_rejects_orphan_mate(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    read1.write_text("@spot/1\nAC\n+\nII\n", encoding="ascii")
    read2.write_text("", encoding="ascii")
    try:
        fastq_pair_stats(read1, read2)
    except ValueError as error:
        assert "count mismatch" in str(error)
    else:
        raise AssertionError("orphan FASTQ mate was accepted")


def test_fastq_pair_stats_requires_mate_labels(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"; read2 = tmp_path / "r2.fastq"
    for path in (read1, read2): path.write_text("@spot\nAC\n+\nII\n", encoding="ascii")
    try:
        fastq_pair_stats(read1, read2)
    except ValueError as error:
        assert "name mismatch" in str(error)
    else:
        raise AssertionError("unlabelled duplicate mates were accepted")
    read1.write_text("@spot 1 length=2\nAC\n+\nII\n", encoding="ascii")
    read2.write_text("@spot 2 length=2\nGT\n+\nII\n", encoding="ascii")
    assert fastq_pair_stats(read1, read2) == {"pair_count": 1, "read_count_per_file": 1, "total_bases": 4}


def test_forbidden_replicate_scan_is_recursive_and_case_insensitive(tmp_path: Path) -> None:
    nested = tmp_path / "nested"; nested.mkdir()
    forbidden = nested / "srr31989027_partial.fastq.gz"; forbidden.write_bytes(b"partial")
    assert forbidden_replicate_files(tmp_path) == [forbidden]


def test_verified_index_prefix_binds_manifest_paths(tmp_path: Path) -> None:
    values = {}
    for index in range(1, 9):
        path = tmp_path / f"strain-b.{index}.ht2"; path.write_bytes(str(index).encode())
        values[path.name] = {"path": path.relative_to(tmp_path).as_posix(), "size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    assert verified_index_prefix(tmp_path, values, "strain-b", "test") == tmp_path / "strain-b"
    values["strain-b.1.ht2"]["path"] = "../outside.ht2"
    try:
        verified_index_prefix(tmp_path, values, "strain-b", "test")
    except ValueError:
        pass
    else:
        raise AssertionError("index path escape was accepted")


def test_splice_support_separates_secondary_and_excludes_supplementary(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"; fasta.write_text(">CP014715.1\n" + "A" * 100 + "\n>CP014716.1\n" + "A" * 100 + "\n", encoding="utf-8")
    sam = "\n".join([
        "q1\t99\tCP014715.1\t1\t60\t5M10N5M\t=\t30\t40\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:1",
        "q1\t147\tCP014715.1\t30\t60\t10M\t=\t1\t-40\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:1",
        "q1\t355\tCP014716.1\t2\t10\t5M8N5M\t=\t30\t40\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:2",
        "q1\t2147\tCP014716.1\t3\t10\t5M7N5M\t=\t30\t40\tAAAAAAAAAA\tFFFFFFFFFF\tNH:i:1",
    ]) + "\n"
    summary = tmp_path / "summary.json"; junctions = tmp_path / "junctions.tsv"
    completed = subprocess.run([sys.executable, str(ROOT / "scripts/summarize_srr31989028_alignments.py"), "--reference-fasta", str(fasta), "--summary-json", str(summary), "--junctions-tsv", str(junctions)], input=sam, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    value = json.loads(summary.read_text(encoding="utf-8"))["junction_summary"]
    assert value == {"distinct_junctions": 2, "with_multi_fragment_support": 1, "with_unique_fragment_support": 1}


def test_acceptance_schema_rejects_capability_overclaim() -> None:
    schema = json.loads((ROOT / "schemas/srr31989028_probe_acceptance.v1.schema.json").read_text(encoding="utf-8"))
    minimal_identity = {"size_bytes": 1, "sha256": "0" * 64}
    value = {
        "schema_version": 1, "acceptance_version": "srr31989028-controlled-probe-v1", "run_id": "srr31989028-probe-" + "0" * 16,
        "probe_id": "prjna1210090-srr31989028-controlled-probe", "probe_acceptance_status": "passed", "execution_status": "complete", "verification_status": "passed", "scientific_acceptance_status": "blocked",
        "capabilities": {"coordinate_coverage": "passed", "splice_support": "passed", "strand_support": "available", "boundary_support": "not-authorized", "strain_specific_support": "unavailable", "biological_replication": "single-replicate-only", "unique_multi_unmapped_classification": "passed", "engineering_screen": "passed-with-explicit-local-exclusions"},
        "strand_support": "unavailable", "single_replicate_only": True, "three_replicate_evidence": False, "formal_boundary_track": False,
        "input_manifest": minimal_identity, "acquisition_evidence": minimal_identity, "source_sra": minimal_identity, "run_manifest": minimal_identity, "artifacts": {str(index): minimal_identity for index in range(16)},
        "verification_evidence": {"independent": minimal_identity, "automated_tests": minimal_identity | {"passed_count": 1}, "repeatability": minimal_identity},
        "protected_evidence": {name: {} for name in ("slice0", "slice0a", "slice0b", "prjna604658", "prjna1210090_metadata")},
        "stop_line": "single SRR31989028 probe complete; remain blocked and do not download other replicates, create a boundary track, generate candidates or thresholds, or enter Slice 1",
    }
    assert list(jsonschema.Draft202012Validator(schema).iter_errors(value))
