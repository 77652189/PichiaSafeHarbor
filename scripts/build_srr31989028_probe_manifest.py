from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import jsonschema


PROTECTED = {
    "slice0": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json",
    "slice0a": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json",
    "slice0b": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json",
    "prjna604658": "local_runs/independent_transcript_sources/prjna604658/qualification_v6_run1/acceptance_manifest.json",
    "prjna1210090_metadata": "local_runs/independent_transcript_sources/prjna1210090/qualification_v4_run1/acceptance_manifest.json",
}
LENGTHS = {"CP014715.1": 2887587, "CP014716.1": 2396528, "CP014717.1": 2249375, "CP014718.1": 1831004, "CP014724.1": 11426, "CP014719.1": 13861, "CP014720.1": 31580, "CP014721.1": 13305, "CP014722.1": 14703, "CP014723.1": 13388}
ENDPOINTS = {"CP014715.1": "5prime_unknown;3prime_partial", "CP014716.1": "complete", "CP014717.1": "5prime_unknown;3prime_partial", "CP014718.1": "5prime_unknown;3prime_partial"}


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def identity(path: Path, root: Path) -> dict:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise ValueError(f"file must be inside repo root: {path}")
    return {"path": resolved.relative_to(root.resolve()).as_posix(), "size_bytes": resolved.stat().st_size, "sha256": digest(resolved, "sha256")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(); root = args.repo_root.resolve(); output = args.output.resolve()
    if not output.is_relative_to(root): raise SystemExit("manifest output must be inside repo root")
    if args.output.exists(): raise SystemExit("probe manifest already exists")
    base = root / "local_runs/controlled_probe/prjna1210090_srr31989028"
    raw = base / "raw"
    acquisition_path = raw / "acquisition_evidence.v1.json"
    if not acquisition_path.is_file(): raise SystemExit("missing acquisition evidence")
    acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
    jsonschema.validate(acquisition, json.loads((root / "schemas/srr31989028_acquisition_evidence.v1.schema.json").read_text(encoding="utf-8")))
    nested_identities = [
        acquisition["source"]["metadata_snapshot"], acquisition["download_receipt"], acquisition["sra"],
        acquisition["validation"]["tool"], acquisition["validation"]["stdout"], acquisition["validation"]["stderr"],
        acquisition["conversion"]["tool"], acquisition["conversion"]["stdout"], acquisition["conversion"]["stderr"],
        acquisition["compression"]["stdout"], acquisition["compression"]["stderr"],
    ]
    for value in nested_identities:
        if identity(root / value["path"], root) != {key: value[key] for key in ("path", "size_bytes", "sha256")}:
            raise SystemExit(f"acquisition artifact identity mismatch: {value['path']}")
    if digest(root / acquisition["sra"]["path"], "md5") != acquisition["sra"]["md5"]:
        raise SystemExit("normalized SRA MD5 mismatch")
    if acquisition.get("acquisition_method") != "ncbi-aws-normalized-sra-to-deterministic-paired-fastq" or acquisition.get("run_accession") != "SRR31989028" or acquisition.get("forbidden_replicates_present") is not False:
        raise SystemExit("acquisition evidence contract mismatch")
    expected = {"read1": "SRR31989028_1.fastq.gz", "read2": "SRR31989028_2.fastq.gz"}
    reads = {}
    for key, name in expected.items():
        path = raw / name
        file_identity = identity(path, root)
        if acquisition.get("reads", {}).get(key) != file_identity:
            raise SystemExit(f"derived FASTQ identity mismatch: {key}")
        reads[key] = file_identity | {"derived_from_sra_sha256": acquisition["sra"]["sha256"], "pair_count": acquisition["conversion"]["statistics"]["pair_count"], "compression": "python-gzip-zlib-level6-mtime0-empty-name"}
    forbidden = [path for path in raw.rglob("*") if path.is_file() and path.name.upper().startswith(("SRR31989016", "SRR31989027"))]
    if forbidden: raise SystemExit("forbidden replicate raw files detected")
    reference_manifest = root / "reference/manifest.v1.json"
    reference = json.loads(reference_manifest.read_text(encoding="utf-8"))["references"]["strain-b"]
    fasta = root / "reference/data/strain-b" / reference["files"]["fasta"]["local_name"]
    gff = root / "reference/data/strain-b" / reference["files"]["annotation"]["local_name"]
    if identity(fasta, root)["sha256"] != reference["files"]["fasta"]["sha256"] or identity(gff, root)["sha256"] != reference["files"]["annotation"]["sha256"]:
        raise SystemExit("reference identity mismatch")
    index_files = sorted((base / "reference/hisat2").glob("strain-b.*.ht2"))
    univec = base / "external_panel/UniVec_Core"; univec_indexes = sorted((base / "external_panel/hisat2").glob("univec.*.ht2"))
    if len(index_files) != 8 or len(univec_indexes) != 8: raise SystemExit("HISAT2 index file set mismatch")
    version_result = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", "/usr/bin/hisat2 --version | head -1; /usr/bin/samtools --version | head -1; dpkg-query -W hisat2 samtools"], capture_output=True, check=True)
    versions = version_result.stdout.decode("utf-8", errors="replace").replace("\x00", "").strip().splitlines()
    windows = []
    for reference_name in ("CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"):
        length = LENGTHS[reference_name]; middle = length // 2
        for position_class, start, end in (("start", 0, 2000), ("middle", middle - 1000, middle + 1000), ("end", length - 2000, length)):
            windows.append({"window_id": f"{reference_name}:{position_class}", "reference": reference_name, "start_0based": start, "end_0based": end, "position_class": position_class, "endpoint_caveat": ENDPOINTS[reference_name] if position_class in {"start", "end"} else "not_terminal"})
    protected = {name: identity(root / relative, root) for name, relative in PROTECTED.items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    windows_path = output.parent / "fixed_windows.v1.json"
    windows_path.write_text(json.dumps({"schema_version": 1, "coordinate_space": "GCA_001746955.1", "sequence_lengths": LENGTHS, "windows": windows}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": 1, "probe_id": "prjna1210090-srr31989028-controlled-probe", "run_accession": "SRR31989028",
        "capability_scope": {"coordinate_coverage": "probe", "splice_support": "probe", "strand_support": "unavailable", "boundary_support": "not-authorized", "strain_specific_support": "identity_needs_sequence_check", "biological_replication": "single-replicate-only"},
        "raw_data_governance": {"local_compute_only": True, "redistribution": False, "publish_raw_reads": False, "archive_terms": "archive-terms-with-submitter-caveat"},
        "acquisition": {"evidence": identity(acquisition_path, root), "sra": acquisition["sra"], "metadata_snapshot": acquisition["source"]["metadata_snapshot"], "validation_status": acquisition["validation"]["status"], "conversion_status": acquisition["conversion"]["status"], "source_url": acquisition["source"]["url"]},
        "reads": reads,
        "reference": {"assembly_accession": "GCA_001746955.1", "manifest": identity(reference_manifest, root), "fasta": identity(fasta, root), "gff": identity(gff, root), "sequence_lengths": LENGTHS},
        "tools": {"wsl_distribution": "Ubuntu-24.04", "versions": versions, "hisat2_index": {path.name: identity(path, root) for path in index_files}, "univec": identity(univec, root), "univec_hisat2_index": {path.name: identity(path, root) for path in univec_indexes}},
        "parameters": {"hisat2": ["-p", "8", "--seed", "42", "--reorder", "--dta", "-k", "5", "--max-intronlen", "10000", "--no-temp-splicesite", "--new-summary"], "samtools_sort": ["--no-PG", "-@", "4", "-m", "1G", "-O", "BAM"], "depth": ["-aa", "-d", "0"], "input_order": ["read1", "read2"], "randomness": "HISAT2 seed 42 with --reorder; non-deterministic mode disabled"},
        "engineering_screening_rules": {"version": "controlled-probe-v1", "global_univec_signal_if_fraction_of_all_pairs_gte": 0.01, "global_univec_signal_if_fraction_of_both_unmapped_pairs_gte": 0.10, "widespread_high_depth_if_fraction_of_nuclear_1kb_bins_gte": 0.05, "interpretation": "probe QC descriptors only; not safe-harbor, expression, or boundary thresholds"},
        "fixed_windows": windows, "fixed_windows_file": identity(windows_path, root), "protected_authority": protected,
        "stop_line": "single SRR31989028 probe only; do not download other replicates, create a boundary track, generate candidates or thresholds, or enter Slice 1",
    }
    schema = json.loads((root / "schemas/srr31989028_probe_manifest.v1.schema.json").read_text(encoding="utf-8")); jsonschema.validate(manifest, schema)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"probe_id": manifest["probe_id"], "read_sha256": {key: value["sha256"] for key, value in reads.items()}, "windows": len(windows)}, indent=2))
    return 0


if __name__ == "__main__": raise SystemExit(main())
