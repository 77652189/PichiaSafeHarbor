from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import shlex
import statistics
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zlib
from collections import Counter, defaultdict
from pathlib import Path

import jsonschema


NUCLEAR = {"CP014715.1": 2887587, "CP014716.1": 2396528, "CP014717.1": 2249375, "CP014718.1": 1831004}
CANONICAL_UNSTRANDED = {"GT-AG", "GC-AG", "AT-AC", "CT-AC", "CT-GC", "GT-AT"}
EXPECTED_ARTIFACTS = {"alignment.bam", "alignment.bam.bai", "hisat2_summary.txt", "hisat2_stderr.txt", "flagstat.json", "samtools_stats.txt", "pair_mapping_summary.json", "splice_junctions.tsv", "coverage_summary.json", "coverage_windows.tsv", "coverage_anomalies.tsv", "univec_all_pairs_summary.txt", "univec_both_unmapped_summary.txt", "engineering_impact.json", "capability_report.json", "capability_report.md"}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
TOOLCHAIN = (
    "scripts/record_srr31989028_download.py", "scripts/prepare_srr31989028_reads.py", "scripts/build_srr31989028_probe_manifest.py", "scripts/run_srr31989028_controlled_probe.py",
    "scripts/summarize_srr31989028_alignments.py", "scripts/summarize_srr31989028_coverage.py",
    "scripts/verify_srr31989028_controlled_probe.py", "scripts/run_srr31989028_probe_workflow.py",
    "schemas/srr31989028_acquisition_evidence.v1.schema.json", "schemas/srr31989028_probe_manifest.v1.schema.json", "schemas/srr31989028_probe_acceptance.v1.schema.json",
    "tests/test_srr31989028_controlled_probe.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def fastq_name_and_mate(header: bytes) -> tuple[bytes, int | None]:
    fields = header.rstrip(b"\r\n").split(); name = fields[0]
    if name.endswith((b"/1", b"/2")): return name[:-2], int(name[-1:])
    if len(fields) > 1 and fields[1] in {b"1", b"2"}: return name, int(fields[1])
    return name, None


def recompute_fastq_stats(read1: Path, read2: Path) -> dict:
    pairs = 0; bases = 0
    with gzip.open(read1, "rb") as left, gzip.open(read2, "rb") as right:
        while True:
            records = []
            for handle in (left, right):
                record = [handle.readline() for _ in range(4)]
                if not record[0]:
                    if any(record): raise SystemExit("truncated compressed FASTQ")
                    records.append(None); continue
                if any(not line for line in record) or not record[0].startswith(b"@") or not record[2].startswith(b"+") or len(record[1].rstrip()) != len(record[3].rstrip()): raise SystemExit("invalid compressed FASTQ record")
                records.append(record)
            if records == [None, None]: break
            if None in records: raise SystemExit("compressed FASTQ pair count mismatch")
            left_name, left_mate = fastq_name_and_mate(records[0][0]); right_name, right_mate = fastq_name_and_mate(records[1][0])
            if left_name != right_name or left_mate != 1 or right_mate != 2: raise SystemExit("compressed FASTQ mate mismatch")
            pairs += 1; bases += len(records[0][1].rstrip()) + len(records[1][1].rstrip())
    return {"pair_count": pairs, "read_count_per_file": pairs, "total_bases": bases}


def verify_metadata_snapshot(path: Path) -> None:
    root = ET.parse(path).getroot(); run = next((node for node in root.iter("RUN") if node.get("accession") == "SRR31989028"), None)
    files = [] if run is None else [node for node in run.iter("SRAFile") if node.get("semantic_name") == "SRA Normalized"]
    if len(files) != 1 or files[0].get("url") != "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR31989028/SRR31989028" or files[0].get("size") != "3931837205" or files[0].get("md5") != "8c8a58880254746890e10c68b402c875": raise SystemExit("independent NCBI metadata identity mismatch")


def identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": sha256(path)}


def wsl_path(path: Path) -> str:
    value = path.resolve(); return f"/mnt/{value.drive[0].lower()}{value.as_posix()[2:]}"


def toolchain(root: Path) -> str:
    digest = hashlib.sha256()
    for relative in TOOLCHAIN:
        path = root / relative; digest.update(relative.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


def runtime_versions() -> list[str]:
    command = "/usr/bin/hisat2 --version | head -1; /usr/bin/samtools --version | head -1; dpkg-query -W hisat2 samtools"
    completed = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], capture_output=True, check=True)
    return completed.stdout.decode("utf-8", errors="replace").replace("\x00", "").strip().splitlines()


def parse_summary(path: Path) -> dict:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*([A-Za-z][A-Za-z >0-9_-]*):\s*([0-9.]+)%?", line)
        if match: values[match.group(1).strip().lower().replace(" ", "_").replace(">", "gt")] = float(match.group(2))
    return values


def verified_index_prefix(root: Path, values: dict, prefix_name: str) -> Path:
    expected = {f"{prefix_name}.{index}.ht2" for index in range(1, 9)}
    if set(values) != expected: raise SystemExit("independent index file set mismatch")
    paths = []
    for name in sorted(values):
        path = (root / values[name]["path"]).resolve()
        if not path.is_relative_to(root) or identity(path) != {"size_bytes": values[name]["size_bytes"], "sha256": values[name]["sha256"]}: raise SystemExit("independent index identity mismatch")
        paths.append(path)
    if len({path.parent for path in paths}) != 1: raise SystemExit("independent index directory mismatch")
    return paths[0].parent / prefix_name


def protected(root: Path, value: dict, label: str) -> dict:
    acceptance_path = (root / value["path"]).resolve()
    if not acceptance_path.is_relative_to(root) or identity(acceptance_path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}: raise SystemExit(f"protected acceptance mismatch: {label}")
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8")); run_identity = acceptance["run_manifest"]; run_path = acceptance_path.parent / "run_manifest.json"
    if identity(run_path) != run_identity: raise SystemExit(f"protected run mismatch: {label}")
    run = json.loads(run_path.read_text(encoding="utf-8")); artifacts = {}
    for name, artifact_identity in sorted(run["artifacts"].items()):
        path = (acceptance_path.parent / name).resolve()
        if not path.is_relative_to(acceptance_path.parent.resolve()) or identity(path) != artifact_identity: raise SystemExit(f"protected artifact mismatch: {label}:{name}")
        artifacts[name] = artifact_identity["sha256"]
    return {"acceptance_sha256": value["sha256"], "run_manifest_sha256": run_identity["sha256"], "artifact_sha256": artifacts}


def nh(row: list[str]) -> int | None:
    for value in row[11:]:
        if value.startswith("NH:i:"): return int(value[5:])
    return None


def load_fasta(path: Path) -> dict[str, str]:
    result = {}; current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]; result[current] = []
        elif current is not None:
            result[current].append(line.strip())
    return {name: "".join(parts).upper() for name, parts in result.items()}


def process(records: list[list[str]], counts: Counter, junctions: dict) -> None:
    if not records: return
    counts["total_query_pairs"] += 1; primary = [row for row in records if not int(row[1]) & (0x100 | 0x800)]
    r1 = next((row for row in primary if int(row[1]) & 0x40), None); r2 = next((row for row in primary if int(row[1]) & 0x80), None)
    if r1 is None or r2 is None: counts["unknown_pair_structure"] += 1
    else:
        mapped = [not int(row[1]) & 0x4 for row in (r1, r2)]
        if mapped == [False, False]: counts["pair_unmapped"] += 1
        elif mapped[0] != mapped[1]: counts["pair_mixed"] += 1
        else:
            counts["pair_both_mapped"] += 1
            if int(r1[1]) & 0x2 and int(r2[1]) & 0x2: counts["pair_proper"] += 1
            else: counts["pair_discordant"] += 1
            tags = [nh(r1), nh(r2)]; secondary = any(int(row[1]) & 0x100 for row in records)
            if tags == [1, 1] and not secondary: counts["pair_unique"] += 1
            elif secondary or any(value is not None and value > 1 for value in tags): counts["pair_multi"] += 1
            else: counts["pair_multiplicity_unknown"] += 1
    unique_seen = set(); multi_seen = set()
    for row in records:
        flag = int(row[1])
        if flag & (0x4 | 0x800) or row[5] == "*": continue
        target = unique_seen if nh(row) == 1 and not (flag & (0x100 | 0x800)) else multi_seen
        ref = int(row[3]) - 1
        for length_text, op in CIGAR_RE.findall(row[5]):
            length = int(length_text)
            if op == "N": target.add((row[2], ref, ref + length)); ref += length
            elif op in {"M", "D", "=", "X"}: ref += length
    for key in unique_seen: junctions[key]["unique"] += 1
    for key in multi_seen: junctions[key]["multi"] += 1


def recompute_pairs(bam: Path) -> tuple[dict, dict, dict]:
    command = f"set -o pipefail; /usr/bin/samtools collate -O -u {shlex.quote(wsl_path(bam))} | /usr/bin/samtools view -h -"
    process_handle = subprocess.Popen(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
    counts = Counter(); junctions = defaultdict(Counter); current = None; group = []
    assert process_handle.stdout is not None
    for line in process_handle.stdout:
        if line.startswith("@"): continue
        row = line.rstrip("\n").split("\t")
        if current is not None and row[0] != current: process(group, counts, junctions); group = []
        current = row[0]; group.append(row)
    process(group, counts, junctions); stderr = process_handle.stderr.read() if process_handle.stderr else ""; code = process_handle.wait()
    if code != 0: raise SystemExit(f"independent pair stream failed: {stderr}")
    summary = {"distinct_junctions": len(junctions), "with_unique_fragment_support": sum(value["unique"] > 0 for value in junctions.values()), "with_multi_fragment_support": sum(value["multi"] > 0 for value in junctions.values())}
    return dict(counts), summary, {key: {"unique": value["unique"], "multi": value["multi"]} for key, value in junctions.items()}


def recompute_coverage(bam: Path, lengths: dict, windows: list[dict]) -> dict:
    command = f"/usr/bin/samtools depth -aa -d 0 {shlex.quote(wsl_path(bam))}"
    handle = subprocess.Popen(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
    result = {name: {"positions": 0, "covered": 0, "sum": 0, "bins": defaultdict(int)} for name in lengths}
    window_result = {value["window_id"]: {"covered": 0, "sum": 0} for value in windows}; windows_by_ref = defaultdict(list)
    for value in windows: windows_by_ref[value["reference"]].append(value)
    assert handle.stdout is not None
    for line in handle.stdout:
        reference, _, depth_text = line.rstrip("\n").split("\t")[:3]
        if reference in result:
            depth = int(depth_text); pos0 = int(line.split("\t", 2)[1]) - 1; result[reference]["positions"] += 1; result[reference]["covered"] += depth > 0; result[reference]["sum"] += depth; result[reference]["bins"][pos0 // 1000] += depth
            for window in windows_by_ref[reference]:
                if window["start_0based"] <= pos0 < window["end_0based"]:
                    window_result[window["window_id"]]["covered"] += depth > 0; window_result[window["window_id"]]["sum"] += depth
    stderr = handle.stderr.read() if handle.stderr else ""; code = handle.wait()
    if code != 0: raise SystemExit(f"independent coverage stream failed: {stderr}")
    for name, length in lengths.items():
        if result[name]["positions"] != length: raise SystemExit(f"independent depth length mismatch: {name}")
    anomaly_rows = []
    for name, length in lengths.items():
        means = [result[name]["bins"].get(index, 0) / min(1000, length - index * 1000) for index in range((length + 999) // 1000)]
        nonzero = [value for value in means if value > 0]; median = statistics.median(nonzero) if nonzero else 0.0; mad = statistics.median(abs(value - median) for value in nonzero) if nonzero else 0.0; threshold = median + max(10 * mad, 5 * median, 10.0)
        for index, value in enumerate(means):
            if value > threshold:
                anomaly_rows.append({"reference": name, "start_0based": index * 1000, "end_0based": min(length, (index + 1) * 1000), "mean_depth": value, "chromosome_median_nonzero_bin_depth": median, "chromosome_mad_nonzero_bin_depth": mad, "classification": "high_depth_outlier"})
        del result[name]["bins"]
    anomaly_rows.sort(key=lambda row: (-row["mean_depth"], row["reference"], row["start_0based"]))
    return {"per_sequence": result, "windows": window_result, "anomaly_rows": anomaly_rows, "high_depth_outlier_bins": len(anomaly_rows), "nuclear_high_depth_outlier_bins": sum(row["reference"] in NUCLEAR for row in anomaly_rows)}


def recompute_univec(root: Path, manifest: dict, bam: Path, read1: Path, read2: Path) -> dict:
    prefix = verified_index_prefix(root, manifest["tools"]["univec_hisat2_index"], "univec")
    with tempfile.TemporaryDirectory(prefix="verify_srr31989028_", dir=bam.parent) as directory:
        temp = Path(directory); u1 = temp / "u1.fastq.gz"; u2 = temp / "u2.fastq.gz"; all_summary = temp / "all.txt"; unmapped_summary = temp / "unmapped.txt"
        extract = f"set -o pipefail; /usr/bin/samtools collate -O -u {shlex.quote(wsl_path(bam))} | /usr/bin/samtools fastq -@ 4 -f 12 -F 2304 -n -c 1 -1 {shlex.quote(wsl_path(u1))} -2 {shlex.quote(wsl_path(u2))} -0 /dev/null -s /dev/null -"
        commands = [extract, f"/usr/bin/hisat2 -p 4 --seed 42 -k 1 --no-spliced-alignment --no-unal --new-summary --summary-file {shlex.quote(wsl_path(all_summary))} -x {shlex.quote(wsl_path(prefix))} -1 {shlex.quote(wsl_path(read1))} -2 {shlex.quote(wsl_path(read2))} -S /dev/null", f"/usr/bin/hisat2 -p 4 --seed 42 -k 1 --no-spliced-alignment --no-unal --new-summary --summary-file {shlex.quote(wsl_path(unmapped_summary))} -x {shlex.quote(wsl_path(prefix))} -1 {shlex.quote(wsl_path(u1))} -2 {shlex.quote(wsl_path(u2))} -S /dev/null"]
        for command in commands:
            completed = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], capture_output=True, text=True, encoding="utf-8", errors="replace", check=False)
            if completed.returncode != 0: raise SystemExit(f"independent UniVec step failed: {completed.stderr}")
        keys = ("aligned_concordantly_1_time", "aligned_concordantly_gt1_times", "aligned_discordantly_1_time")
        all_values = parse_summary(all_summary); unmapped_values = parse_summary(unmapped_summary)
        return {"all": int(sum(all_values.get(key, 0) for key in keys)), "both_unmapped": int(sum(unmapped_values.get(key, 0) for key in keys)), "all_total_pairs": int(all_values.get("total_pairs", -1)), "both_unmapped_total_pairs": int(unmapped_values.get("total_pairs", -1))}


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("--run-dir", type=Path, required=True); parser.add_argument("--manifest", type=Path, required=True); parser.add_argument("--repo-root", type=Path, default=Path(".")); parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(); root = args.repo_root.resolve(); manifest_path = args.manifest.resolve(); run_dir = args.run_dir.resolve(); output = args.output.resolve()
    if not manifest_path.is_relative_to(root) or not manifest_path.is_file(): raise SystemExit("manifest must be a file inside repo root")
    if not run_dir.is_relative_to(root) or not run_dir.is_dir(): raise SystemExit("run directory must be inside repo root")
    if not output.is_relative_to(run_dir): raise SystemExit("verification output must be inside run directory")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")); jsonschema.validate(manifest, json.loads((root / "schemas/srr31989028_probe_manifest.v1.schema.json").read_text(encoding="utf-8"))); run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8")); workflow_hash = toolchain(root)
    if run.get("input_manifest_sha256") != sha256(manifest_path) or run.get("toolchain_sha256") != workflow_hash: raise SystemExit("manifest or toolchain mismatch")
    if manifest["reference"]["assembly_accession"] != "GCA_001746955.1" or manifest["capability_scope"]["strand_support"] != "unavailable": raise SystemExit("coordinate or strand scope mismatch")
    raw = root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw"
    if any(path.is_file() and path.name.upper().startswith(("SRR31989016", "SRR31989027")) for path in raw.rglob("*")): raise SystemExit("forbidden replicate raw files detected")
    if manifest["tools"].get("wsl_distribution") != "Ubuntu-24.04" or runtime_versions() != manifest["tools"].get("versions"): raise SystemExit("independent runtime tool version mismatch")
    read_paths = []
    for read in manifest["reads"].values():
        path = (root / read["path"]).resolve()
        if not path.is_relative_to(root) or identity(path) != {"size_bytes": read["size_bytes"], "sha256": read["sha256"]}: raise SystemExit("FASTQ identity mismatch")
        read_paths.append(path)
    acquisition_value = manifest["acquisition"]["evidence"]
    acquisition_path = (root / acquisition_value["path"]).resolve()
    if not acquisition_path.is_relative_to(root) or identity(acquisition_path) != {"size_bytes": acquisition_value["size_bytes"], "sha256": acquisition_value["sha256"]}:
        raise SystemExit("acquisition evidence identity mismatch")
    acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
    jsonschema.validate(acquisition, json.loads((root / "schemas/srr31989028_acquisition_evidence.v1.schema.json").read_text(encoding="utf-8")))
    expected_acquired_reads = {key: {field: manifest["reads"][key][field] for field in ("path", "size_bytes", "sha256")} for key in ("read1", "read2")}
    if acquisition.get("reads") != expected_acquired_reads or acquisition.get("validation", {}).get("status") != "passed" or acquisition.get("conversion", {}).get("status") != "passed":
        raise SystemExit("independent acquisition cross-check mismatch")
    nested_identities = [
        acquisition["source"]["metadata_snapshot"], acquisition["download_receipt"], acquisition["sra"],
        acquisition["validation"]["tool"], acquisition["validation"]["stdout"], acquisition["validation"]["stderr"],
        acquisition["conversion"]["tool"], acquisition["conversion"]["stdout"], acquisition["conversion"]["stderr"],
        acquisition["compression"]["stdout"], acquisition["compression"]["stderr"],
    ]
    for value in nested_identities:
        path = (root / value["path"]).resolve()
        if not path.is_relative_to(root) or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}:
            raise SystemExit(f"independent acquisition artifact mismatch: {value['path']}")
    source_sra = (root / acquisition["sra"]["path"]).resolve()
    if md5(source_sra) != acquisition["sra"]["md5"]:
        raise SystemExit("independent normalized SRA MD5 mismatch")
    metadata_path = (root / acquisition["source"]["metadata_snapshot"]["path"]).resolve(); verify_metadata_snapshot(metadata_path)
    download_receipt = json.loads((root / acquisition["download_receipt"]["path"]).read_text(encoding="utf-8"))
    actual_sra = {"size_bytes": source_sra.stat().st_size, "md5": md5(source_sra), "sha256": sha256(source_sra)}
    if download_receipt.get("download_status") != "complete" or download_receipt.get("source_url") != acquisition["source"]["url"] or download_receipt.get("declared_size_bytes") != 3931837205 or download_receipt.get("declared_md5") != acquisition["source"]["declared_md5"] or download_receipt.get("actual_sra") != actual_sra or download_receipt.get("tool_version") != "aria2 version 1.37.0": raise SystemExit("independent download receipt mismatch")
    aria2_value = download_receipt.get("tool", {}); aria2_path = (root / aria2_value.get("path", "")).resolve()
    if not aria2_path.is_relative_to(root) or identity(aria2_path) != {"size_bytes": aria2_value.get("size_bytes"), "sha256": aria2_value.get("sha256")}: raise SystemExit("independent aria2 tool identity mismatch")
    validation_argv = acquisition["validation"]["argv"]; conversion_argv = acquisition["conversion"]["argv"]
    expected_validator = str((root / acquisition["validation"]["tool"]["path"]).resolve())
    compression = acquisition["compression"]
    if validation_argv != [expected_validator, str(source_sra)] or conversion_argv[-1] != str(source_sra) or conversion_argv[1:8] != ["--split-files", "-e", "8", "--seq-defline", "@$ac.$si $ri length=$rl", "--qual-defline", "+$ac.$si $ri length=$rl"] or "-t" not in conversion_argv or "-O" not in conversion_argv or compression["tool"] != "python-standard-library:gzip.GzipFile" or compression["version"] != f"Python {platform.python_version()}; zlib {zlib.ZLIB_VERSION}" or compression["parameters"] != ["compresslevel=6", "mtime=0", "filename="]: raise SystemExit("independent acquisition argv mismatch")
    fastq_stats = recompute_fastq_stats(read_paths[0], read_paths[1])
    if fastq_stats != {"pair_count": 21957417, "read_count_per_file": 21957417, "total_bases": 6587225100} or fastq_stats != acquisition["conversion"]["statistics"]: raise SystemExit("independent FASTQ statistics mismatch")
    verified_reference = {}
    for key in ("manifest", "fasta", "gff"):
        value = manifest["reference"][key]; path = (root / value["path"]).resolve()
        if not path.is_relative_to(root) or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}: raise SystemExit(f"reference identity mismatch: {key}")
        verified_reference[key] = path
    if set(run.get("artifacts", {})) != EXPECTED_ARTIFACTS: raise SystemExit("artifact allowlist mismatch")
    artifacts = {}
    for name, value in run["artifacts"].items():
        path = (run_dir / name).resolve()
        if not path.is_relative_to(run_dir): raise SystemExit(f"artifact path escape: {name}")
        if identity(path) != value: raise SystemExit(f"artifact mismatch: {name}")
        artifacts[name] = value["sha256"]
    bam = run_dir / "alignment.bam"
    header = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", f"/usr/bin/samtools view -H {shlex.quote(wsl_path(bam))}"], capture_output=True, text=True, encoding="utf-8", errors="replace", check=True).stdout.replace("\x00", "")
    header_lengths = {match.group(1): int(match.group(2)) for match in re.finditer(r"@SQ\tSN:([^\t]+)\tLN:(\d+)", header)}
    if any(header_lengths.get(name) != length for name, length in manifest["reference"]["sequence_lengths"].items()): raise SystemExit("BAM coordinate header mismatch")
    pair_counts, junction_summary, junction_support = recompute_pairs(bam); pair_report = json.loads((run_dir / "pair_mapping_summary.json").read_text(encoding="utf-8")); keys = {"total_query_pairs", "pair_unmapped", "pair_mixed", "pair_both_mapped", "pair_unique", "pair_multi", "pair_multiplicity_unknown", "unknown_pair_structure"}
    if any(pair_counts.get(key, 0) != pair_report["counts"].get(key, 0) for key in keys) or junction_summary != pair_report["junction_summary"]: raise SystemExit("independent pair or splice summary mismatch")
    sequences = load_fasta(verified_reference["fasta"]); expected_junction_rows = {}
    for (reference, start, end), support in junction_support.items():
        sequence = sequences.get(reference, ""); motif = f"{sequence[start:start+2]}-{sequence[end-2:end]}" if 0 <= start < end <= len(sequence) else "unavailable"
        expected_junction_rows[(reference, start, end)] = {"motif_unstranded": motif, "canonical_unstranded": str(motif in CANONICAL_UNSTRANDED).lower(), "unique_fragment_support": support["unique"], "multi_fragment_support": support["multi"], "strand_support": "unavailable"}
    with (run_dir / "splice_junctions.tsv").open(encoding="utf-8", newline="") as handle:
        reported_junction_rows = {(row["reference"], int(row["intron_start_0based"]), int(row["intron_end_0based"])): {"motif_unstranded": row["motif_unstranded"], "canonical_unstranded": row["canonical_unstranded"], "unique_fragment_support": int(row["unique_fragment_support"]), "multi_fragment_support": int(row["multi_fragment_support"]), "strand_support": row["strand_support"]} for row in csv.DictReader(handle, delimiter="\t")}
    if reported_junction_rows != expected_junction_rows: raise SystemExit("independent splice junction table mismatch")
    fixed_path = (root / manifest["fixed_windows_file"]["path"]).resolve()
    if not fixed_path.is_relative_to(root) or identity(fixed_path) != {"size_bytes": manifest["fixed_windows_file"]["size_bytes"], "sha256": manifest["fixed_windows_file"]["sha256"]}: raise SystemExit("fixed windows identity mismatch")
    fixed = json.loads(fixed_path.read_text(encoding="utf-8"))
    if fixed.get("sequence_lengths") != manifest["reference"]["sequence_lengths"] or fixed.get("windows") != manifest["fixed_windows"]: raise SystemExit("fixed windows inline/file mismatch")
    coverage = recompute_coverage(bam, manifest["reference"]["sequence_lengths"], fixed["windows"]); coverage_report = json.loads((run_dir / "coverage_summary.json").read_text(encoding="utf-8"))
    for name, value in coverage["per_sequence"].items():
        reported = coverage_report["per_sequence"][name]
        if value["covered"] != reported["covered_bases"] or abs(value["sum"] / value["positions"] - reported["mean_depth"]) > 1e-12: raise SystemExit(f"independent coverage mismatch: {name}")
    if coverage["high_depth_outlier_bins"] != coverage_report["high_depth_outlier_bins"] or coverage["nuclear_high_depth_outlier_bins"] != coverage_report["nuclear_high_depth_outlier_bins"]: raise SystemExit("independent coverage outlier mismatch")
    with (run_dir / "coverage_windows.tsv").open(encoding="utf-8", newline="") as handle:
        reported_windows = {row["window_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(reported_windows) != set(coverage["windows"]): raise SystemExit("fixed window set mismatch")
    for window in fixed["windows"]:
        expected = coverage["windows"][window["window_id"]]; reported = reported_windows[window["window_id"]]; width = window["end_0based"] - window["start_0based"]
        if int(reported["covered_bases"]) != expected["covered"] or abs(float(reported["mean_depth"]) - expected["sum"] / width) > 1e-12: raise SystemExit(f"fixed window mismatch: {window['window_id']}")
    protected_evidence = {name: protected(root, value, name) for name, value in manifest["protected_authority"].items()}
    if protected_evidence != run.get("protected_evidence"): raise SystemExit("protected evidence mismatch")
    capability = json.loads((run_dir / "capability_report.json").read_text(encoding="utf-8")); engineering = json.loads((run_dir / "engineering_impact.json").read_text(encoding="utf-8"))
    if capability.get("strand_direction") != "unavailable" or capability.get("single_replicate_only") is not True or capability.get("three_replicate_evidence") is not False or capability.get("formal_boundary_track") is not False: raise SystemExit("capability overclaim")
    if engineering.get("unresolved_engineering_effects") is None: raise SystemExit("engineering limitations missing")
    univec = recompute_univec(root, manifest, bam, read_paths[0], read_paths[1]); metrics = engineering["metrics"]
    if univec["all"] != metrics["univec_aligned_pairs_all_reads"] or univec["both_unmapped"] != metrics["univec_aligned_pairs_both_unmapped"]: raise SystemExit("independent UniVec mismatch")
    total = pair_counts.get("total_query_pairs", 0); both_unmapped = pair_counts.get("pair_unmapped", 0); nuclear_bins = sum(value["total_1kb_bins"] for value in coverage_report["per_sequence"].values() if value["class"] == "nuclear_chromosome")
    if univec["all_total_pairs"] != total or univec["both_unmapped_total_pairs"] != both_unmapped: raise SystemExit("independent UniVec input pair count mismatch")
    with (run_dir / "coverage_anomalies.tsv").open(encoding="utf-8", newline="") as handle:
        reported_anomalies = list(csv.DictReader(handle, delimiter="\t"))
    expected_anomalies = [{key: str(row[key]) for key in ("reference", "start_0based", "end_0based", "mean_depth", "chromosome_median_nonzero_bin_depth", "chromosome_mad_nonzero_bin_depth", "classification")} for row in coverage["anomaly_rows"]]
    if reported_anomalies != expected_anomalies: raise SystemExit("independent coverage anomaly table mismatch")
    localized = [row for row in expected_anomalies if row["reference"] in NUCLEAR]
    fractions = {"all": univec["all"] / total if total else 0, "unmapped": univec["both_unmapped"] / both_unmapped if both_unmapped else 0, "outlier": len(localized) / nuclear_bins if nuclear_bins else 0}
    rules = manifest["engineering_screening_rules"]; global_signal = fractions["all"] >= rules["global_univec_signal_if_fraction_of_all_pairs_gte"] or fractions["unmapped"] >= rules["global_univec_signal_if_fraction_of_both_unmapped_pairs_gte"]; widespread = fractions["outlier"] >= rules["widespread_high_depth_if_fraction_of_nuclear_1kb_bins_gte"]
    if engineering.get("global_univec_signal") != global_signal or engineering.get("widespread_high_depth") != widespread or engineering.get("localized_suspect_regions") != localized: raise SystemExit("independent engineering gate mismatch")
    expected_capabilities = {"coordinate_coverage": "passed" if all(coverage["per_sequence"][name]["covered"] > 0 for name in NUCLEAR) else "failed", "splice_support": "passed" if junction_summary["with_unique_fragment_support"] > 0 else "failed", "strand_support": "unavailable", "boundary_support": "not-authorized", "strain_specific_support": "unavailable", "biological_replication": "single-replicate-only", "unique_multi_unmapped_classification": "passed" if pair_report["classification_complete"] and pair_counts.get("pair_multiplicity_unknown", 0) == 0 else "partial", "engineering_screen": "passed-with-explicit-local-exclusions" if not (global_signal or widespread) else "failed"}
    if capability.get("capabilities") != expected_capabilities: raise SystemExit("independent capability derivation mismatch")
    hisat_values = parse_summary(run_dir / "hisat2_summary.txt")
    if int(hisat_values.get("total_pairs", -1)) != total: raise SystemExit("HISAT2/BAM total pair mismatch")
    result = {"schema_version": 1, "verification_type": "srr31989028_controlled_probe", "status": "passed", "run_id": run["run_id"], "input_manifest_sha256": run["input_manifest_sha256"], "toolchain_sha256": workflow_hash, "acquisition_evidence_sha256": acquisition_value["sha256"], "source_sra_sha256": acquisition["sra"]["sha256"], "recomputed_fastq_statistics": fastq_stats, "verified_artifacts": artifacts, "recomputed_pair_counts": {key: pair_counts.get(key, 0) for key in sorted(keys)}, "recomputed_junction_summary": junction_summary, "recomputed_coverage": coverage, "recomputed_univec": univec, "protected_evidence": protected_evidence, "strand_support": "unavailable", "single_replicate_only": True, "scientific_acceptance_status": "blocked"}
    output.parent.mkdir(parents=True, exist_ok=True); output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"); print(json.dumps({"status": "passed", "run_id": run["run_id"], "pairs": pair_counts.get("total_query_pairs"), "junctions": junction_summary["distinct_junctions"]}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
