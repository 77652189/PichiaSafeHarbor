from __future__ import annotations

import shutil
import zipfile
import io
from copy import deepcopy
from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.reference import install_reference_archive, verify_reference_files
from pichia_safe_harbor import reference as reference_module


def _archive_and_reference(
    tmp_path: Path, fixture_paths, fixture_reference
) -> tuple[Path, dict]:
    fasta_path, gff_path = fixture_paths
    reference = deepcopy(fixture_reference)
    reference["key"] = "fixture"
    reference["files"]["fasta"].update(
        {"local_name": "fixture.fna", "archive_path": "bundle/fixture.fna"}
    )
    reference["files"]["annotation"].update(
        {"local_name": "fixture.gff3", "archive_path": "bundle/fixture.gff3"}
    )
    archive = tmp_path / "bundle.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.write(fasta_path, "bundle/fixture.fna")
        bundle.write(gff_path, "bundle/fixture.gff3")
    return archive, reference


def _seed_valid_bundle(data_dir: Path, reference: dict, fixture_paths) -> Path:
    target = data_dir / reference["key"]
    target.mkdir(parents=True)
    fasta_path, gff_path = fixture_paths
    shutil.copyfile(fasta_path, target / reference["files"]["fasta"]["local_name"])
    shutil.copyfile(gff_path, target / reference["files"]["annotation"]["local_name"])
    (target / "old-only-sentinel.txt").write_text("old bundle\n", encoding="utf-8")
    verify_reference_files(reference, data_dir)
    return target


@pytest.mark.parametrize(
    "failure_injection",
    ["after_first_file", "after_validation", "after_backup", "after_swap"],
)
def test_atomic_bundle_failure_preserves_previous_valid_bundle(
    tmp_path: Path, fixture_paths, fixture_reference, failure_injection: str
) -> None:
    archive, reference = _archive_and_reference(tmp_path, fixture_paths, fixture_reference)
    data_dir = tmp_path / "data"
    target = _seed_valid_bundle(data_dir, reference, fixture_paths)
    before = {
        path.name: path.read_bytes() for path in target.iterdir() if path.is_file()
    }
    with pytest.raises(ContractError, match="injected"):
        install_reference_archive(
            reference,
            data_dir,
            archive,
            failure_injection=failure_injection,
        )
    after = {path.name: path.read_bytes() for path in target.iterdir() if path.is_file()}
    assert after == before
    verify_reference_files(reference, data_dir)


def test_atomic_bundle_success_replaces_whole_directory(
    tmp_path: Path, fixture_paths, fixture_reference
) -> None:
    archive, reference = _archive_and_reference(tmp_path, fixture_paths, fixture_reference)
    data_dir = tmp_path / "data"
    target = _seed_valid_bundle(data_dir, reference, fixture_paths)
    install_reference_archive(reference, data_dir, archive)
    assert not (target / "old-only-sentinel.txt").exists()
    assert sorted(path.name for path in target.iterdir()) == ["fixture.fna", "fixture.gff3"]
    verify_reference_files(reference, data_dir)


def test_download_retries_from_a_clean_file(monkeypatch, tmp_path: Path) -> None:
    calls = 0

    def fake_urlopen(_request, timeout):
        nonlocal calls
        calls += 1
        if calls < 3:
            raise OSError("injected transport failure")
        return io.BytesIO(b"complete archive")

    monkeypatch.setattr(reference_module.urllib.request, "urlopen", fake_urlopen)
    archive = tmp_path / "archive.zip"
    archive.write_bytes(b"partial")
    reference_module._download_archive("https://example.invalid/archive", archive)
    assert calls == 3
    assert archive.read_bytes() == b"complete archive"
