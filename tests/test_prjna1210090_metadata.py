from __future__ import annotations

import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from scripts.run_prjna1210090_metadata import CANDIDATE_ID, GATES, qualify, verified_file


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "transcript_source_candidates/prjna1210090/manifest.v1.json"


def test_prjna1210090_metadata_qualification_stops_before_raw_reads(tmp_path: Path) -> None:
    output = tmp_path / "run"
    qualify(MANIFEST, output, ROOT)
    result = json.loads((output / "metadata_qualification.json").read_text(encoding="utf-8"))
    assert result["candidate_id"] == CANDIDATE_ID
    assert result["gates"] == GATES
    assert result["gates"]["raw_file_license"] == "unresolved"
    assert result["gates"]["strand_specificity"] == "unavailable"
    assert result["gates"]["biological_replication"] == "passed"
    assert result["gates"]["unmodified_strain-b_identity"] == "unresolved"
    assert result["raw_reads"]["status"] == "not-acquired"
    assert result["eligible_as_independent_transcript_source"] is False
    assert result["scientific_acceptance_status"] == "blocked"


def test_prjna1210090_acquisition_evidence_is_a_hard_identity_gate(tmp_path: Path) -> None:
    evidence = tmp_path / "acquisition.json"
    evidence.write_text("original", encoding="utf-8")
    value = {"path": "acquisition.json", "size_bytes": evidence.stat().st_size, "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest()}
    assert verified_file(tmp_path, value, "acquisition") == evidence
    evidence.write_text("tampered", encoding="utf-8")
    with pytest.raises(ValueError, match="identity mismatch"):
        verified_file(tmp_path, value, "acquisition")


def test_prjna1210090_acquisition_path_cannot_escape_repo(tmp_path: Path) -> None:
    outside = tmp_path.parent / "outside.json"
    outside.write_text("outside", encoding="utf-8")
    value = {"path": "../outside.json", "size_bytes": outside.stat().st_size, "sha256": hashlib.sha256(outside.read_bytes()).hexdigest()}
    with pytest.raises(ValueError, match="identity mismatch"):
        verified_file(tmp_path, value, "acquisition")


def test_prjna1210090_manifest_schema() -> None:
    schema = json.loads((ROOT / "schemas/prjna1210090_metadata_manifest.v1.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(json.loads(MANIFEST.read_text(encoding="utf-8")), schema)


def test_prjna1210090_selected_runs_are_three_distinct_biosamples(tmp_path: Path) -> None:
    output = tmp_path / "run"
    qualify(MANIFEST, output, ROOT)
    result = json.loads((output / "metadata_qualification.json").read_text(encoding="utf-8"))
    assert {item["replicate"] for item in result["records"]} == {"1", "2", "3"}
    assert len({item["sample"] for item in result["records"]}) == 3


def test_prjna1210090_policy_is_not_promoted_to_file_license(tmp_path: Path) -> None:
    output = tmp_path / "run"
    qualify(MANIFEST, output, ROOT)
    result = json.loads((output / "metadata_qualification.json").read_text(encoding="utf-8"))
    assert result["license_evidence"]["status"] == "unresolved"
    assert "do not grant" in result["license_evidence"]["interpretation"]


def test_prjna1210090_reference_manifest_identity_is_required(tmp_path: Path) -> None:
    value = json.loads(MANIFEST.read_text(encoding="utf-8"))
    value["files"]["reference_manifest"]["sha256"] = "0" * 64
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="reference_manifest"):
        qualify(manifest, tmp_path / "run", ROOT)


def test_prjna1210090_protected_authority_cannot_be_replaced(tmp_path: Path) -> None:
    value = json.loads(MANIFEST.read_text(encoding="utf-8"))
    value["protected_authority"] = {}
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="protected authority"):
        qualify(manifest, tmp_path / "run", ROOT)
