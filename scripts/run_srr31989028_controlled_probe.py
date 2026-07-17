from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

import jsonschema


EXPECTED_ARTIFACTS = {
    "alignment.bam", "alignment.bam.bai", "hisat2_summary.txt", "hisat2_stderr.txt", "flagstat.json", "samtools_stats.txt",
    "pair_mapping_summary.json", "splice_junctions.tsv", "coverage_summary.json", "coverage_windows.tsv", "coverage_anomalies.tsv",
    "univec_all_pairs_summary.txt", "univec_both_unmapped_summary.txt", "engineering_impact.json", "capability_report.json", "capability_report.md",
}
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


def identity(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": sha256(path)}


def wsl_path(path: Path) -> str:
    value = path.resolve()
    return f"/mnt/{value.drive[0].lower()}{value.as_posix()[2:]}"


def toolchain_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    for relative in TOOLCHAIN:
        path = root / relative
        if not path.is_file(): raise ValueError(f"toolchain file missing: {relative}")
        digest.update(relative.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


def runtime_versions() -> list[str]:
    command = "/usr/bin/hisat2 --version | head -1; /usr/bin/samtools --version | head -1; dpkg-query -W hisat2 samtools"
    completed = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], capture_output=True, check=True)
    return completed.stdout.decode("utf-8", errors="replace").replace("\x00", "").strip().splitlines()


def verified_file(root: Path, value: dict, label: str) -> Path:
    path = (root.resolve() / value["path"]).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file() or identity(path) != {"size_bytes": value["size_bytes"], "sha256": value["sha256"]}:
        raise ValueError(f"input identity mismatch: {label}")
    return path


def verified_index_prefix(root: Path, values: dict, prefix_name: str, label: str) -> Path:
    expected = {f"{prefix_name}.{index}.ht2" for index in range(1, 9)}
    if set(values) != expected:
        raise ValueError(f"index file set mismatch: {label}")
    paths = [verified_file(root, values[name], f"{label}:{name}") for name in sorted(values)]
    parents = {path.parent for path in paths}
    if len(parents) != 1:
        raise ValueError(f"index directory mismatch: {label}")
    return paths[0].parent / prefix_name


def verify_protected(root: Path, value: dict, label: str) -> dict:
    acceptance_path = verified_file(root, value, label); acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    run_identity = acceptance["run_manifest"]; run_path = acceptance_path.parent / "run_manifest.json"
    if identity(run_path) != run_identity: raise ValueError(f"protected run mismatch: {label}")
    run = json.loads(run_path.read_text(encoding="utf-8")); artifacts = {}
    for name, artifact_identity in sorted(run["artifacts"].items()):
        path = (acceptance_path.parent / name).resolve()
        if not path.is_relative_to(acceptance_path.parent.resolve()) or identity(path) != artifact_identity: raise ValueError(f"protected artifact mismatch: {label}:{name}")
        artifacts[name] = artifact_identity["sha256"]
    return {"acceptance_sha256": value["sha256"], "run_manifest_sha256": run_identity["sha256"], "artifact_sha256": artifacts}


def run_wsl(command: str, cwd: Path) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", command], cwd=cwd, capture_output=True, check=False)
    ended = datetime.now(timezone.utc).isoformat()
    evidence = {"command": command, "started_at_utc": started, "completed_at_utc": ended, "exit_code": completed.returncode, "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(), "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest()}
    if completed.returncode != 0: raise RuntimeError(f"WSL command failed ({completed.returncode}): {command}\n{completed.stdout.decode('utf-8', errors='replace')}\n{completed.stderr.decode('utf-8', errors='replace')}")
    return evidence


def parse_summary(path: Path) -> dict:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*([A-Za-z][A-Za-z >0-9_-]*):\s*([0-9.]+)%?", line)
        if match: values[match.group(1).strip().lower().replace(" ", "_").replace(">", "gt")] = float(match.group(2))
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True); parser.add_argument("--output-dir", type=Path, required=True); parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args(); root = args.repo_root.resolve(); output_dir = args.output_dir.resolve()
    manifest_path = args.manifest.resolve()
    if not manifest_path.is_relative_to(root) or not manifest_path.is_file(): raise SystemExit("manifest must be a file inside repo root")
    if not output_dir.is_relative_to(root): raise SystemExit("output directory must be inside repo root")
    if output_dir.exists(): raise SystemExit(f"output exists: {output_dir}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    jsonschema.validate(manifest, json.loads((root / "schemas/srr31989028_probe_manifest.v1.schema.json").read_text(encoding="utf-8")))
    if manifest.get("probe_id") != "prjna1210090-srr31989028-controlled-probe" or manifest.get("run_accession") != "SRR31989028" or manifest["capability_scope"]["strand_support"] != "unavailable": raise SystemExit("probe manifest contract mismatch")
    raw = root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw"
    if any(path.is_file() and path.name.upper().startswith(("SRR31989016", "SRR31989027")) for path in raw.rglob("*")): raise SystemExit("forbidden replicate raw files detected")
    if manifest["tools"].get("wsl_distribution") != "Ubuntu-24.04" or runtime_versions() != manifest["tools"].get("versions"): raise SystemExit("runtime tool version mismatch")
    read1 = verified_file(root, manifest["reads"]["read1"], "read1"); read2 = verified_file(root, manifest["reads"]["read2"], "read2")
    acquisition_path = verified_file(root, manifest["acquisition"]["evidence"], "acquisition evidence")
    acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
    expected_acquired_reads = {key: {field: manifest["reads"][key][field] for field in ("path", "size_bytes", "sha256")} for key in ("read1", "read2")}
    if acquisition.get("reads") != expected_acquired_reads or acquisition.get("validation", {}).get("status") != "passed" or acquisition.get("conversion", {}).get("status") != "passed":
        raise SystemExit("acquisition/read cross-check mismatch")
    source_sra = verified_file(root, manifest["acquisition"]["sra"], "normalized SRA")
    verified_file(root, manifest["acquisition"]["metadata_snapshot"], "NCBI metadata snapshot")
    if md5(source_sra) != manifest["acquisition"]["sra"]["md5"]:
        raise SystemExit("normalized SRA MD5 mismatch")
    fasta = verified_file(root, manifest["reference"]["fasta"], "reference fasta"); verified_file(root, manifest["reference"]["gff"], "reference gff"); verified_file(root, manifest["reference"]["manifest"], "reference manifest")
    index_prefix_path = verified_index_prefix(root, manifest["tools"]["hisat2_index"], "strain-b", "hisat2 index")
    univec = verified_file(root, manifest["tools"]["univec"], "UniVec")
    univec_prefix_path = verified_index_prefix(root, manifest["tools"]["univec_hisat2_index"], "univec", "UniVec index")
    protected = {name: verify_protected(root, value, name) for name, value in manifest["protected_authority"].items()}
    workflow_hash = toolchain_sha256(root); material = json.dumps({"manifest_sha256": sha256(manifest_path), "toolchain_sha256": workflow_hash}, sort_keys=True).encode(); run_id = "srr31989028-probe-" + hashlib.sha256(material).hexdigest()[:16]
    output_dir.parent.mkdir(parents=True, exist_ok=True); temp = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent)); steps = []
    try:
        paths = {name: wsl_path(temp / name) for name in EXPECTED_ARTIFACTS}
        index_prefix = wsl_path(index_prefix_path)
        align = "set -o pipefail; " + " ".join([
            "/usr/bin/hisat2", *map(shlex.quote, manifest["parameters"]["hisat2"]), "--rg-id", "SRR31989028", "--rg", "SM:SRR31989028",
            "--summary-file", shlex.quote(paths["hisat2_summary.txt"]), "-x", shlex.quote(index_prefix), "-1", shlex.quote(wsl_path(read1)), "-2", shlex.quote(wsl_path(read2)),
            "2>" + shlex.quote(paths["hisat2_stderr.txt"]), "|", "/usr/bin/awk", shlex.quote('substr($0,1,3)!="@PG"'), "|", "/usr/bin/samtools", "sort", *map(shlex.quote, manifest["parameters"]["samtools_sort"]), "-o", shlex.quote(paths["alignment.bam"]), "-",
        ])
        steps.append(run_wsl(align, root))
        steps.append(run_wsl(f"/usr/bin/samtools index -@ 4 {shlex.quote(paths['alignment.bam'])} {shlex.quote(paths['alignment.bam.bai'])}", root))
        steps.append(run_wsl(f"/usr/bin/samtools flagstat -O json {shlex.quote(paths['alignment.bam'])} > {shlex.quote(paths['flagstat.json'])} && /usr/bin/samtools stats - < {shlex.quote(paths['alignment.bam'])} > {shlex.quote(paths['samtools_stats.txt'])}", root))
        summary_script = wsl_path(root / "scripts/summarize_srr31989028_alignments.py")
        steps.append(run_wsl(f"set -o pipefail; /usr/bin/samtools collate -O -u {shlex.quote(paths['alignment.bam'])} | /usr/bin/samtools view -h - | /usr/bin/python3 {shlex.quote(summary_script)} --reference-fasta {shlex.quote(wsl_path(fasta))} --summary-json {shlex.quote(paths['pair_mapping_summary.json'])} --junctions-tsv {shlex.quote(paths['splice_junctions.tsv'])}", root))
        windows_json = verified_file(root, manifest["fixed_windows_file"], "fixed windows")
        coverage_script = wsl_path(root / "scripts/summarize_srr31989028_coverage.py")
        steps.append(run_wsl(f"set -o pipefail; /usr/bin/samtools depth -aa -d 0 {shlex.quote(paths['alignment.bam'])} | /usr/bin/python3 {shlex.quote(coverage_script)} --windows-json {shlex.quote(wsl_path(windows_json))} --summary-json {shlex.quote(paths['coverage_summary.json'])} --windows-tsv {shlex.quote(paths['coverage_windows.tsv'])} --anomalies-tsv {shlex.quote(paths['coverage_anomalies.tsv'])}", root))
        unmapped1 = temp / "unmapped_1.fastq.gz"; unmapped2 = temp / "unmapped_2.fastq.gz"
        steps.append(run_wsl(f"set -o pipefail; /usr/bin/samtools collate -O -u {shlex.quote(paths['alignment.bam'])} | /usr/bin/samtools fastq -@ 4 -f 12 -F 2304 -n -c 1 -1 {shlex.quote(wsl_path(unmapped1))} -2 {shlex.quote(wsl_path(unmapped2))} -0 /dev/null -s /dev/null -", root))
        univec_prefix = wsl_path(univec_prefix_path)
        steps.append(run_wsl(f"/usr/bin/hisat2 -p 4 --seed 42 -k 1 --no-spliced-alignment --no-unal --new-summary --summary-file {shlex.quote(paths['univec_all_pairs_summary.txt'])} -x {shlex.quote(univec_prefix)} -1 {shlex.quote(wsl_path(read1))} -2 {shlex.quote(wsl_path(read2))} -S /dev/null", root))
        steps.append(run_wsl(f"/usr/bin/hisat2 -p 4 --seed 42 -k 1 --no-spliced-alignment --no-unal --new-summary --summary-file {shlex.quote(paths['univec_both_unmapped_summary.txt'])} -x {shlex.quote(univec_prefix)} -1 {shlex.quote(wsl_path(unmapped1))} -2 {shlex.quote(wsl_path(unmapped2))} -S /dev/null", root))
        unmapped1.unlink(missing_ok=True); unmapped2.unlink(missing_ok=True)
        pair = json.loads((temp / "pair_mapping_summary.json").read_text(encoding="utf-8")); coverage = json.loads((temp / "coverage_summary.json").read_text(encoding="utf-8")); univec_all = parse_summary(temp / "univec_all_pairs_summary.txt"); univec_unmapped = parse_summary(temp / "univec_both_unmapped_summary.txt")
        counts = pair["counts"]; total = counts["total_query_pairs"]; both_unmapped = counts.get("pair_unmapped", 0)
        aligned_keys = ("aligned_concordantly_1_time", "aligned_concordantly_gt1_times", "aligned_discordantly_1_time")
        univec_all_aligned = int(sum(univec_all.get(key, 0) for key in aligned_keys)); univec_unmapped_aligned = int(sum(univec_unmapped.get(key, 0) for key in aligned_keys))
        if int(univec_all.get("total_pairs", -1)) != total or int(univec_unmapped.get("total_pairs", -1)) != both_unmapped: raise RuntimeError("UniVec input pair count mismatch")
        nuclear_bins = sum(value["total_1kb_bins"] for value in coverage["per_sequence"].values() if value["class"] == "nuclear_chromosome")
        rules = manifest["engineering_screening_rules"]; univec_all_fraction = univec_all_aligned / total if total else 0; univec_unmapped_fraction = univec_unmapped_aligned / both_unmapped if both_unmapped else 0; anomaly_fraction = coverage["nuclear_high_depth_outlier_bins"] / nuclear_bins if nuclear_bins else 0
        with (temp / "coverage_anomalies.tsv").open(encoding="utf-8", newline="") as handle:
            localized_regions = [row for row in csv.DictReader(handle, delimiter="\t") if row["reference"] in {"CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"}]
        engineering = {
            "schema_version": 1, "screening_scope": "UniVec/adaptor-like sequence, non-reference burden, discordant/supplementary/soft-clipped pairs, and robust 1kb coverage outliers",
            "global_univec_signal": univec_all_fraction >= rules["global_univec_signal_if_fraction_of_all_pairs_gte"] or univec_unmapped_fraction >= rules["global_univec_signal_if_fraction_of_both_unmapped_pairs_gte"],
            "widespread_high_depth": anomaly_fraction >= rules["widespread_high_depth_if_fraction_of_nuclear_1kb_bins_gte"],
            "metrics": {"univec_aligned_pairs_all_reads": univec_all_aligned, "univec_aligned_pairs_both_unmapped": univec_unmapped_aligned, "univec_fraction_all_pairs": univec_all_fraction, "univec_fraction_both_unmapped_pairs": univec_unmapped_fraction, "high_depth_outlier_fraction_nuclear_bins": anomaly_fraction, "pair_softclip_ge20": counts.get("pair_softclip_ge20", 0), "pair_indel_ge20": counts.get("pair_indel_ge20", 0), "pair_with_supplementary": counts.get("pair_with_supplementary", 0), "pair_discordant": counts.get("pair_discordant", 0), "pair_non_nuclear_reference": counts.get("pair_non_nuclear_reference", 0)},
            "localized_suspect_regions": localized_regions,
            "unresolved_engineering_effects": ["UniVec is not a comprehensive catalog of engineered genes or integration constructs", "RNA-seq high-depth and soft-clipped loci are expression/breakpoint candidates, not proof of genomic integration", "single RNA replicate cannot exclude low-expression or silent engineered sequence", "flagged coverage outliers are unavailable for strain-specific boundary inference"],
            "capability_conclusion": "no global signal" if not (univec_all_fraction >= rules["global_univec_signal_if_fraction_of_all_pairs_gte"] or univec_unmapped_fraction >= rules["global_univec_signal_if_fraction_of_both_unmapped_pairs_gte"] or anomaly_fraction >= rules["widespread_high_depth_if_fraction_of_nuclear_1kb_bins_gte"]) else "global or widespread signal detected",
        }
        (temp / "engineering_impact.json").write_text(json.dumps(engineering, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        capabilities = {
            "coordinate_coverage": "passed" if coverage["nuclear_chromosomes_nonzero"] else "failed", "splice_support": "passed" if pair["junction_summary"]["with_unique_fragment_support"] > 0 else "failed",
            "strand_support": "unavailable", "boundary_support": "not-authorized", "strain_specific_support": "unavailable", "biological_replication": "single-replicate-only",
            "unique_multi_unmapped_classification": "passed" if pair["classification_complete"] and counts.get("pair_multiplicity_unknown", 0) == 0 else "partial",
            "engineering_screen": "passed-with-explicit-local-exclusions" if engineering["capability_conclusion"] == "no global signal" else "failed",
        }
        capability = {"schema_version": 1, "run_id": run_id, "probe_id": manifest["probe_id"], "scientific_acceptance_status": "blocked", "capabilities": capabilities, "strand_direction": "unavailable", "single_replicate_only": True, "three_replicate_evidence": False, "formal_boundary_track": False, "stop_line": manifest["stop_line"]}
        (temp / "capability_report.json").write_text(json.dumps(capability, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (temp / "capability_report.md").write_text("# SRR31989028 controlled probe\n\n- Results are limited to one WT-labelled replicate and Strain-B GCA_001746955.1 coordinates.\n- Strand direction remains unavailable.\n- Coverage, splice, unique, multi, unmapped, and engineering-impact screens are reported separately.\n- This is not three-replicate evidence, a formal boundary track, a safe-harbor candidate run, or a threshold decision.\n", encoding="utf-8")
        artifacts = {name: identity(temp / name) for name in sorted(EXPECTED_ARTIFACTS)}
        run_manifest = {"schema_version": 1, "run_id": run_id, "probe_id": manifest["probe_id"], "execution_status": "complete", "verification_status": "not_run", "scientific_acceptance_status": "blocked", "input_manifest_sha256": sha256(manifest_path), "toolchain_sha256": workflow_hash, "parameters": manifest["parameters"], "artifacts": artifacts, "protected_evidence": protected}
        (temp / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        receipt = {"schema_version": 1, "evidence_type": "controlled_probe_execution", "run_id": run_id, "probe_id": manifest["probe_id"], "invocation_id": str(uuid.uuid4()), "input_manifest_sha256": run_manifest["input_manifest_sha256"], "toolchain_sha256": workflow_hash, "steps": steps, "completed_at_utc": datetime.now(timezone.utc).isoformat()}
        receipt_path = temp / "verification/execution_receipt.json"; receipt_path.parent.mkdir(parents=True); receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temp, output_dir)
        print(json.dumps({"run_id": run_id, "capabilities": capabilities, "bam_size_bytes": artifacts["alignment.bam"]["size_bytes"]}, indent=2)); return 0
    except Exception:
        shutil.rmtree(temp, ignore_errors=True); raise


if __name__ == "__main__": raise SystemExit(main())
