from __future__ import annotations

import json
import http.client
import os
import shutil
import tempfile
import urllib.request
import urllib.error
import uuid
import zipfile
from pathlib import Path
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file


def load_reference_manifest(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("schema_version") != 1:
        raise ContractError("unsupported reference manifest schema_version")
    if not manifest.get("references"):
        raise ContractError("reference manifest contains no references")
    return manifest


def reference_by_key(manifest: dict[str, Any], key: str) -> dict[str, Any]:
    try:
        return manifest["references"][key]
    except KeyError as exc:
        raise ContractError(f"unknown reference key: {key}") from exc


def verify_reference_files(reference: dict[str, Any], data_dir: Path) -> dict[str, Path]:
    resolved = {
        role: data_dir / reference["key"] / reference["files"][role]["local_name"]
        for role in ("fasta", "annotation")
    }
    return verify_reference_paths(reference, resolved["fasta"], resolved["annotation"])


def verify_reference_paths(
    reference: dict[str, Any], fasta_path: Path, annotation_path: Path
) -> dict[str, Path]:
    resolved = {"fasta": fasta_path, "annotation": annotation_path}
    for role, path in resolved.items():
        spec = reference["files"][role]
        if not path.is_file():
            raise ContractError(f"missing {role} file: {path}")
        size = path.stat().st_size
        if size != spec["size_bytes"]:
            raise ContractError(
                f"{role} size mismatch for {path}: expected {spec['size_bytes']}, got {size}"
            )
        digest = sha256_file(path)
        if digest != spec["sha256"]:
            raise ContractError(
                f"{role} sha256 mismatch for {path}: expected {spec['sha256']}, got {digest}"
            )
    return resolved


def fetch_reference(reference: dict[str, Any], data_dir: Path) -> dict[str, Path]:
    """Download an NCBI Datasets archive and install only manifest-declared files."""
    with tempfile.TemporaryDirectory(prefix="pichia-safe-harbor-") as tmp_name:
        tmp = Path(tmp_name)
        archive = tmp / "ncbi_dataset.zip"
        _download_archive(reference["download"]["url"], archive)
        return install_reference_archive(reference, data_dir, archive)


def _download_archive(url: str, archive: Path, attempts: int = 3) -> None:
    last_error: Exception | None = None
    for _attempt in range(1, attempts + 1):
        archive.unlink(missing_ok=True)
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": "PichiaSafeHarbor/0.1",
                "Connection": "close",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response, archive.open(
                "wb"
            ) as out:
                shutil.copyfileobj(response, out)
            return
        except (OSError, urllib.error.URLError, http.client.HTTPException) as exc:
            last_error = exc
    raise ContractError(f"reference download failed after {attempts} attempts: {last_error}")


def install_reference_archive(
    reference: dict[str, Any],
    data_dir: Path,
    archive: Path,
    *,
    failure_injection: str | None = None,
) -> dict[str, Path]:
    """Install a complete verified bundle without exposing a mixed-version directory."""
    data_dir.mkdir(parents=True, exist_ok=True)
    data_root = data_dir.resolve()
    target_dir = (data_dir / reference["key"]).resolve()
    if target_dir.parent != data_root:
        raise ContractError(f"reference target escapes data directory: {target_dir}")

    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{reference['key']}.staging-", dir=data_dir)
    ).resolve()
    backup_dir = (data_dir / f".{reference['key']}.backup-{uuid.uuid4().hex}").resolve()
    backup_created = False
    try:
        with zipfile.ZipFile(archive) as bundle:
            names = set(bundle.namelist())
            staged_paths: dict[str, Path] = {}
            for index, role in enumerate(("fasta", "annotation"), 1):
                spec = reference["files"][role]
                member = spec["archive_path"]
                if member not in names:
                    raise ContractError(f"archive is missing declared {role}: {member}")
                destination = staging_dir / spec["local_name"]
                with bundle.open(member) as source, destination.open("wb") as output:
                    shutil.copyfileobj(source, output)
                staged_paths[role] = destination
                if failure_injection == "after_first_file" and index == 1:
                    raise ContractError("injected reference installation failure after first file")

        verify_reference_paths(
            reference, staged_paths["fasta"], staged_paths["annotation"]
        )
        if failure_injection == "after_validation":
            raise ContractError("injected reference installation failure after validation")

        if target_dir.exists():
            os.replace(target_dir, backup_dir)
            backup_created = True
        try:
            if failure_injection == "after_backup":
                raise ContractError("injected reference installation failure after backup")
            os.replace(staging_dir, target_dir)
            if failure_injection == "after_swap":
                raise ContractError("injected reference installation failure after swap")
            installed = verify_reference_files(reference, data_dir)
        except Exception:
            failed_dir = (
                data_dir / f".{reference['key']}.failed-{uuid.uuid4().hex}"
            ).resolve()
            if target_dir.exists():
                os.replace(target_dir, failed_dir)
            if backup_created:
                os.replace(backup_dir, target_dir)
                backup_created = False
            if failed_dir.exists():
                shutil.rmtree(failed_dir, ignore_errors=True)
            raise
        if backup_created:
            shutil.rmtree(backup_dir, ignore_errors=True)
            backup_created = False
        return installed
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir, ignore_errors=True)
        if backup_created and backup_dir.exists() and not target_dir.exists():
            os.replace(backup_dir, target_dir)
