from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import jsonschema


ACCESSION = "SRR31989028"
SOURCE_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR31989028/SRR31989028"
SRA_SIZE = 3_931_837_205
SRA_MD5 = "8c8a58880254746890e10c68b402c875"
EXPECTED_PAIRS = 21_957_417
EXPECTED_BASES = 6_587_225_100


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
    return {
        "path": resolved.relative_to(root.resolve()).as_posix(),
        "size_bytes": resolved.stat().st_size,
        "sha256": digest(resolved, "sha256"),
    }


def fingerprint(path: Path) -> dict:
    return {"size_bytes": path.stat().st_size, "sha256": digest(path, "sha256")}


def metadata_sra_identity(path: Path) -> dict:
    root = ET.parse(path).getroot()
    run = next((node for node in root.iter("RUN") if node.get("accession") == ACCESSION), None)
    if run is None:
        raise ValueError("SRR31989028 missing from NCBI metadata snapshot")
    files = [node for node in run.iter("SRAFile") if node.get("semantic_name") == "SRA Normalized"]
    if len(files) != 1:
        raise ValueError("normalized SRA metadata entry mismatch")
    value = files[0]
    result = {"url": value.get("url"), "size_bytes": int(value.get("size", "-1")), "md5": value.get("md5")}
    if result != {"url": SOURCE_URL, "size_bytes": SRA_SIZE, "md5": SRA_MD5}:
        raise ValueError("normalized SRA metadata identity mismatch")
    return result


def fastq_pair_stats(read1: Path, read2: Path) -> dict:
    pairs = 0
    bases = 0
    with read1.open("rb") as left, read2.open("rb") as right:
        while True:
            records = []
            for handle in (left, right):
                record = [handle.readline() for _ in range(4)]
                if not record[0]:
                    if any(record):
                        raise ValueError("truncated FASTQ record")
                    records.append(None)
                    continue
                if any(not line for line in record) or not record[0].startswith(b"@") or not record[2].startswith(b"+"):
                    raise ValueError("invalid FASTQ record")
                if len(record[1].rstrip(b"\r\n")) != len(record[3].rstrip(b"\r\n")):
                    raise ValueError("FASTQ sequence/quality length mismatch")
                records.append(record)
            if records == [None, None]:
                break
            if None in records:
                raise ValueError("paired FASTQ count mismatch")
            left_name, left_mate = fastq_name_and_mate(records[0][0])
            right_name, right_mate = fastq_name_and_mate(records[1][0])
            if left_name != right_name or left_mate != 1 or right_mate != 2:
                raise ValueError("paired FASTQ name mismatch")
            pairs += 1
            bases += len(records[0][1].rstrip(b"\r\n")) + len(records[1][1].rstrip(b"\r\n"))
    return {"pair_count": pairs, "read_count_per_file": pairs, "total_bases": bases}


def fastq_name_and_mate(header: bytes) -> tuple[bytes, int | None]:
    fields = header.rstrip(b"\r\n").split()
    name = fields[0]
    if name.endswith((b"/1", b"/2")):
        return name[:-2], int(name[-1:])
    if len(fields) > 1 and fields[1] in {b"1", b"2"}:
        return name, int(fields[1])
    return name, None


def forbidden_replicate_files(raw: Path) -> list[Path]:
    prefixes = ("SRR31989016", "SRR31989027")
    return [path for path in raw.rglob("*") if path.is_file() and path.name.upper().startswith(prefixes)]


def run(command: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> dict:
    started = datetime.now(timezone.utc).isoformat()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=False)
    result = {
        "argv": command,
        "cwd": str(cwd),
        "started_at_utc": started,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": completed.returncode,
        "stdout": stdout_path,
        "stderr": stderr_path,
    }
    if completed.returncode != 0:
        stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace")[-4000:]
        raise RuntimeError(f"command failed ({completed.returncode}): {command}\n{stderr_text}")
    return result


def version(executable: Path) -> str:
    completed = subprocess.run([str(executable), "--version"], capture_output=True, check=True)
    output = completed.stdout.decode("utf-8", errors="replace").replace("\x00", "").strip()
    match = re.search(r"(\d+\.\d+\.\d+)\s*$", output)
    if match is None:
        raise ValueError(f"unrecognized tool version output: {output}")
    return match.group(1)


def to_wsl(path: Path) -> str:
    value = path.resolve()
    return f"/mnt/{value.drive[0].lower()}{value.as_posix()[2:]}"


def deterministic_gzip(source: Path) -> Path:
    output = source.with_suffix(source.suffix + ".gz")
    with source.open("rb") as input_handle, output.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, compresslevel=6, mtime=0) as compressed:
            shutil.copyfileobj(input_handle, compressed, length=8 * 1024 * 1024)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--sra", type=Path, required=True)
    parser.add_argument("--metadata-snapshot", type=Path, required=True)
    parser.add_argument("--download-receipt", type=Path, required=True)
    parser.add_argument("--toolkit-bin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    raw = (root / "local_runs/controlled_probe/prjna1210090_srr31989028/raw").resolve()
    sra = args.sra.resolve()
    output = args.output.resolve()
    metadata = args.metadata_snapshot.resolve()
    download_receipt = args.download_receipt.resolve()
    if not raw.is_relative_to(root) or not sra.is_relative_to(raw) or not output.is_relative_to(raw):
        raise SystemExit("acquisition paths must stay inside the controlled raw directory")
    if output.exists() or forbidden_replicate_files(raw):
        raise SystemExit("output exists or forbidden replicate raw files detected")
    if sra.stat().st_size != SRA_SIZE or digest(sra, "md5") != SRA_MD5:
        raise SystemExit("normalized SRA identity mismatch")
    metadata_sra_identity(metadata)
    receipt = json.loads(download_receipt.read_text(encoding="utf-8"))
    actual_sra = {"size_bytes": sra.stat().st_size, "md5": digest(sra, "md5"), "sha256": digest(sra, "sha256")}
    if receipt.get("source_url") != SOURCE_URL or receipt.get("declared_size_bytes") != SRA_SIZE or receipt.get("declared_md5") != SRA_MD5 or receipt.get("download_status") != "complete" or receipt.get("actual_sra") != actual_sra or receipt.get("tool_version") != "aria2 version 1.37.0":
        raise SystemExit("download receipt mismatch")
    toolkit = args.toolkit_bin.resolve()
    validator = toolkit / "vdb-validate.exe"
    fasterq = toolkit / "fasterq-dump.exe"
    if version(validator) != "3.4.1" or version(fasterq) != "3.4.1":
        raise SystemExit("SRA Toolkit version mismatch")
    logs = raw / "acquisition_logs"
    for stale in raw.glob("srr31989028-convert-*"):
        if stale.is_dir(): shutil.rmtree(stale)
    if logs.exists(): shutil.rmtree(logs)
    stage = Path(tempfile.mkdtemp(prefix="srr31989028-convert-", dir=raw))
    logs.mkdir()
    try:
        validation = run([str(validator), str(sra)], root, logs / "vdb_validate.stdout.txt", logs / "vdb_validate.stderr.txt")
        temp = stage / "tmp"
        converted = stage / "out"
        temp.mkdir()
        converted.mkdir()
        conversion = run([str(fasterq), "--split-files", "-e", "8", "--seq-defline", "@$ac.$si $ri length=$rl", "--qual-defline", "+$ac.$si $ri length=$rl", "-t", str(temp), "-O", str(converted), str(sra)], root, logs / "fasterq_dump.stdout.txt", logs / "fasterq_dump.stderr.txt")
        read1 = converted / f"{ACCESSION}_1.fastq"
        read2 = converted / f"{ACCESSION}_2.fastq"
        if {path.name for path in converted.iterdir() if path.is_file()} != {read1.name, read2.name}:
            raise RuntimeError("fasterq-dump output file set mismatch")
        stats = fastq_pair_stats(read1, read2)
        if stats != {"pair_count": EXPECTED_PAIRS, "read_count_per_file": EXPECTED_PAIRS, "total_bases": EXPECTED_BASES}:
            raise RuntimeError(f"FASTQ statistics mismatch: {stats}")
        uncompressed = {"read1": fingerprint(read1), "read2": fingerprint(read2)}
        compression_stdout = logs / "gzip.stdout.txt"; compression_stderr = logs / "gzip.stderr.txt"
        compression_stdout.write_bytes(b""); compression_stderr.write_bytes(b"")
        compression_started = datetime.now(timezone.utc).isoformat()
        with ThreadPoolExecutor(max_workers=2) as pool:
            compressed_reads = list(pool.map(deterministic_gzip, (read1, read2)))
        compression = {"argv": [sys.executable, "gzip.GzipFile", "compresslevel=6", "mtime=0", "filename="], "cwd": str(root), "started_at_utc": compression_started, "completed_at_utc": datetime.now(timezone.utc).isoformat(), "exit_code": 0, "stdout": compression_stdout, "stderr": compression_stderr}
        final1 = raw / f"{ACCESSION}_1.fastq.gz"
        final2 = raw / f"{ACCESSION}_2.fastq.gz"
        if final1.exists() or final2.exists():
            raise RuntimeError("final FASTQ already exists")
        os.replace(compressed_reads[0], final1)
        os.replace(compressed_reads[1], final2)
        evidence = {
            "schema_version": 1,
            "acquisition_method": "ncbi-aws-normalized-sra-to-deterministic-paired-fastq",
            "run_accession": ACCESSION,
            "source": {"url": SOURCE_URL, "declared_size_bytes": SRA_SIZE, "declared_md5": SRA_MD5, "metadata_snapshot": identity(metadata, root)},
            "download_receipt": identity(download_receipt, root),
            "sra": identity(sra, root) | {"md5": digest(sra, "md5"), "format": "SRA Normalized", "lite": False},
            "validation": {key: value for key, value in validation.items() if key not in {"stdout", "stderr"}} | {"tool": identity(validator, root), "version": "3.4.1", "stdout": identity(validation["stdout"], root), "stderr": identity(validation["stderr"], root), "status": "passed"},
            "conversion": {key: value for key, value in conversion.items() if key not in {"stdout", "stderr"}} | {"tool": identity(fasterq, root), "version": "3.4.1", "stdout": identity(conversion["stdout"], root), "stderr": identity(conversion["stderr"], root), "threads": 8, "uncompressed_reads": uncompressed, "statistics": stats, "status": "passed"},
            "compression": {key: value for key, value in compression.items() if key not in {"stdout", "stderr"}} | {"tool": "python-standard-library:gzip.GzipFile", "version": f"Python {platform.python_version()}; zlib {zlib.ZLIB_VERSION}", "parameters": ["compresslevel=6", "mtime=0", "filename="], "stdout": identity(compression["stdout"], root), "stderr": identity(compression["stderr"], root), "status": "passed"},
            "reads": {"read1": identity(final1, root), "read2": identity(final2, root)},
            "forbidden_replicates_present": False,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        schema = json.loads((root / "schemas/srr31989028_acquisition_evidence.v1.schema.json").read_text(encoding="utf-8"))
        jsonschema.validate(evidence, schema)
        output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"status": "passed", "pairs": stats["pair_count"], "reads": evidence["reads"]}, indent=2))
        return 0
    except Exception:
        for path in (raw / f"{ACCESSION}_1.fastq.gz", raw / f"{ACCESSION}_2.fastq.gz"):
            path.unlink(missing_ok=True)
        shutil.rmtree(logs, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(stage, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
