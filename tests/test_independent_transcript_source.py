from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pytest

import pichia_safe_harbor.independent_transcript_source as independent_source_module

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.independent_transcript_source import (
    REPEATABLE_FILES,
    create_independent_source_acceptance,
    qualify_independent_transcript_source,
    toolchain_sha256,
)
from pichia_safe_harbor.pipeline import _implementation_hash


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _identity(path: Path) -> dict:
    return {"path": path.as_posix(), "size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    archive = tmp_path / "archive.sra"
    archive.write_bytes(b"candidate-archive")
    archive_md5 = hashlib.md5(archive.read_bytes()).hexdigest()
    xml = tmp_path / "sra.xml"
    _write(
        xml,
        f"""<EXPERIMENT_PACKAGE_SET>
<EXPERIMENT_PACKAGE><EXPERIMENT accession="SRX7669278"><STUDY_REF accession="SRP247014"/><DESIGN><DESIGN_DESCRIPTION>WT expressing eGFP</DESIGN_DESCRIPTION><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN><PLATFORM><ILLUMINA><INSTRUMENT_MODEL>HiSeq X Ten</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT><SAMPLE><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMN13978086</EXTERNAL_ID></IDENTIFIERS><DESCRIPTION>WT</DESCRIPTION><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>strain</TAG><VALUE>Strain-B</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>collection_date</TAG><VALUE>2019-05-11</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>sample_type</TAG><VALUE>cell culture</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE><RUN_SET><RUN accession="SRR11011828" total_spots="5748446" total_bases="850102478"><SRAFiles><SRAFile semantic_name="SRA Normalized" size="{archive.stat().st_size}" md5="{archive_md5}"/></SRAFiles></RUN></RUN_SET></EXPERIMENT_PACKAGE>
<EXPERIMENT_PACKAGE><EXPERIMENT accession="SRX7669279"><DESIGN><DESIGN_DESCRIPTION>WT/eGFP with Bcy1 overexpression</DESIGN_DESCRIPTION><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN></EXPERIMENT><SAMPLE><DESCRIPTION>Bcy1 overexpression</DESCRIPTION><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>strain</TAG><VALUE>Strain-B-Bcy1</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE><RUN_SET><RUN accession="SRR11011827" total_spots="1" total_bases="1"/></RUN_SET></EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>""",
    )
    ena = tmp_path / "ena.json"
    _write(ena, '[{"run_accession":"SRR11011828","study_accession":"PRJNA604658","sample_accession":"SAMN13978086","experiment_accession":"SRX7669278"}]\n')
    pubmed = tmp_path / "pubmed.xml"
    _write(pubmed, '<PubmedArticleSet><PubmedData><ArticleIdList><ArticleId IdType="pubmed">32737719</ArticleId><ArticleId IdType="doi">10.1007/s10529-020-02977-z</ArticleId></ArticleIdList></PubmedData></PubmedArticleSet>\n')
    bioproject = tmp_path / "bioproject.xml"
    _write(bioproject, '<RecordSet><DocumentSummary><ArchiveID accession="PRJNA604658"/></DocumentSummary></RecordSet>\n')
    crossref = tmp_path / "crossref.json"
    _write(crossref, json.dumps({"message": {"DOI": "10.1007/s10529-020-02977-z", "published": {"date-parts": [[2020, 7, 31]]}, "license": [{"URL": "https://www.springer.com/tdm"}]}}))
    files = {
        "sra_normalized": archive,
        "sra_experiment_packages": xml,
        "ena_read_run": ena,
        "pubmed_publication": pubmed,
        "bioproject": bioproject,
        "crossref_publication": crossref,
    }
    protected = {}
    for name in ("slice0", "slice0a", "slice0b"):
        protected_dir = tmp_path / name
        artifact = protected_dir / "evidence.json"
        _write(artifact, json.dumps({"fixture": name}))
        artifact_identity = _identity(artifact)
        protected_run = protected_dir / "run_manifest.json"
        _write(protected_run, json.dumps({"artifacts": {"evidence.json": {key: artifact_identity[key] for key in ("size_bytes", "sha256")}}}))
        protected_run_identity = _identity(protected_run)
        acceptance = protected_dir / "acceptance_manifest.json"
        _write(acceptance, json.dumps({"run_manifest": {key: protected_run_identity[key] for key in ("size_bytes", "sha256")}}))
        protected[name] = _identity(acceptance) | {"path": acceptance.relative_to(tmp_path).as_posix()}
    independent_source_module.EXPECTED_PROTECTED_AUTHORITY = protected
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "candidate_id": "prjna604658-srr11011828",
        "accessions": {"bioproject": "PRJNA604658", "study": "SRP247014", "experiment": "SRX7669278", "sample": "SAMN13978086", "run": "SRR11011828", "publication_doi": "10.1007/s10529-020-02977-z", "pmid": "32737719"},
        "protected_authority": protected,
        "files": {name: _identity(path) for name, path in files.items()},
    }), encoding="utf-8")
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"references": {"strain-b": {"assembly_accession": "GCA_001746955.1", "annotation_accession": "GCA_001746955.1 GenBank annotation release 2016-09-21"}}}), encoding="utf-8")
    return manifest, reference, tmp_path


def _receipt(run_dir: Path, invocation_id: str, manifest: Path, reference: Path, repo_root: Path) -> None:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    path = run_dir / "verification/execution_receipt.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    steps = [
        {"command": [sys.executable, str(root / "scripts/run_independent_source_qualification.py"), "--manifest", str(manifest.resolve()), "--reference-manifest", str(reference.resolve()), "--output-dir", str(run_dir.resolve()), "--repo-root", str(repo_root.resolve())], "started_at_utc": "2026-07-15T00:00:00+00:00", "completed_at_utc": "2026-07-15T00:00:01+00:00", "exit_code": 0, "stdout_sha256": "a"},
        {"command": [sys.executable, str(root / "scripts/run_test_evidence.py"), "--output", str((run_dir / "verification/test_evidence.json").resolve())], "started_at_utc": "2026-07-15T00:00:01+00:00", "completed_at_utc": "2026-07-15T00:00:02+00:00", "exit_code": 0, "stdout_sha256": "b"},
        {"command": [sys.executable, str(root / "scripts/verify_independent_source.py"), "--run-dir", str(run_dir.resolve()), "--manifest", str(manifest.resolve()), "--reference-manifest", str(reference.resolve()), "--repo-root", str(repo_root.resolve()), "--output", str((run_dir / "verification/independent_check.json").resolve())], "started_at_utc": "2026-07-15T00:00:02+00:00", "completed_at_utc": "2026-07-15T00:00:03+00:00", "exit_code": 0, "stdout_sha256": "c"},
    ]
    path.write_text(json.dumps({
        "evidence_type": "independent_execution_receipt",
        "candidate_id": run["candidate_id"],
        "run_id": run["run_id"],
        "invocation_id": invocation_id,
        "command": steps[0]["command"],
        "manifest_sha256": run["manifest_sha256"],
        "reference_manifest_sha256": run["reference_manifest_sha256"],
        "implementation_sha256": run["implementation_sha256"],
        "toolchain_sha256": run["toolchain_sha256"],
        "workflow_evidence_type": "orchestrated_independent_source_execution",
        "workflow_status": "passed",
        "workflow_cwd": str(root),
        "coordinator_id": "fixture-coordinator",
        "run_slot": "run-a" if run_dir.name == "run" else "run-b",
        "workflow_steps": steps,
        "workflow_outputs": {relative: {"size_bytes": (run_dir / relative).stat().st_size, "sha256": hashlib.sha256((run_dir / relative).read_bytes()).hexdigest()} for relative in REPEATABLE_FILES},
    }), encoding="utf-8")


def _verification(run_dir: Path, peer: Path, tmp_path: Path, manifest: Path, reference: Path, root: Path) -> tuple[Path, Path, Path]:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    qualify_independent_transcript_source(manifest, reference, peer, root)
    candidate_manifest = json.loads(manifest.read_text(encoding="utf-8"))
    protected = {}
    for name, identity in candidate_manifest["protected_authority"].items():
        acceptance_path = root / identity["path"]
        acceptance_value = json.loads(acceptance_path.read_text(encoding="utf-8"))
        protected_run = json.loads((acceptance_path.parent / "run_manifest.json").read_text(encoding="utf-8"))
        protected[name] = {
            "acceptance_sha256": identity["sha256"],
            "run_manifest_sha256": acceptance_value["run_manifest"]["sha256"],
            "artifact_sha256": {key: value["sha256"] for key, value in protected_run["artifacts"].items()},
        }
    independent = run_dir / "verification/independent_check.json"
    independent.parent.mkdir(parents=True, exist_ok=True)
    qualification = json.loads((run_dir / "source_qualification.json").read_text(encoding="utf-8"))
    independent.write_text(json.dumps({"verification_type": "independent_transcript_source_qualification", "status": "passed", "run_id": run["run_id"], "manifest_sha256": run["manifest_sha256"], "reference_manifest_sha256": run["reference_manifest_sha256"], "toolchain_sha256": run["toolchain_sha256"], "verified_artifacts": {name: identity["sha256"] for name, identity in run["artifacts"].items()}, "verified_source_files": {name: identity["sha256"] for name, identity in run["source_files"].items()}, "recomputed_gates": qualification["gates"], "eligible_as_independent_transcript_source": False, "scientific_acceptance_status": "blocked", "protected_evidence": protected}), encoding="utf-8")
    tests = run_dir / "verification/test_evidence.json"
    tests.write_text(json.dumps({"schema_version": 1, "evidence_type": "automated_tests", "status": "passed", "command": ["python", "-m", "pytest", "-q"], "passed_count": 70, "implementation_sha256": _implementation_hash(), "toolchain_sha256": toolchain_sha256(Path(__file__).resolve().parents[1])}), encoding="utf-8")
    (peer / "verification").mkdir(parents=True, exist_ok=True)
    shutil.copy2(independent, peer / "verification/independent_check.json")
    shutil.copy2(tests, peer / "verification/test_evidence.json")
    _receipt(run_dir, "invocation-run", manifest, reference, root)
    _receipt(peer, "invocation-peer", manifest, reference, root)
    repeat = tmp_path / "repeatability.json"
    hashes = {relative: hashlib.sha256((run_dir / relative).read_bytes()).hexdigest() for relative in REPEATABLE_FILES}
    receipts = {}
    for current in (run_dir, peer):
        path = current / "verification/execution_receipt.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        receipts[current.name] = {"invocation_id": value["invocation_id"], "identity": {"size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}}
    repeat.write_text(json.dumps({"verification_type": "independent_transcript_source_repeatability", "status": "passed", "run_id": run["run_id"], "toolchain_sha256": run["toolchain_sha256"], "all_files_identical": True, "run_directories": [run_dir.name, peer.name], "file_sha256": hashes, "execution_receipts": receipts}), encoding="utf-8")
    return independent, tests, repeat


def test_candidate_qualification_preserves_blocking_gates(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    result = json.loads((run_dir / "source_qualification.json").read_text(encoding="utf-8"))
    assert result["gates"]["annotation_generation_independence"] == "passed"
    assert result["gates"]["raw_file_license"] == "unresolved"
    assert result["gates"]["strand_specificity"] == "unavailable"
    assert result["gates"]["biological_replication"] == "not-passed"
    assert result["eligible_as_independent_transcript_source"] is False


def test_candidate_qualification_rejects_tampered_source_file(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    (tmp_path / "archive.sra").write_bytes(b"tampered")
    with pytest.raises(ContractError):
        qualify_independent_transcript_source(manifest, reference, tmp_path / "run", root)


def test_candidate_qualification_rejects_wrong_primary_reference(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    value = json.loads(reference.read_text(encoding="utf-8"))
    value["references"]["strain-b"]["assembly_accession"] = "GCA_WRONG"
    reference.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ContractError):
        qualify_independent_transcript_source(manifest, reference, tmp_path / "run", root)


def test_candidate_qualification_rejects_replaced_protected_authority(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    value = json.loads(manifest.read_text(encoding="utf-8"))
    value["protected_authority"] = {}
    manifest.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ContractError):
        qualify_independent_transcript_source(manifest, reference, tmp_path / "run", root)


def test_candidate_acceptance_requires_blocked_qualification_and_real_peer(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    accepted = create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)
    assert accepted["verification_status"] == "passed"
    assert accepted["scientific_acceptance_status"] == "blocked"
    assert accepted["blocking_gates"]["unmodified_strain-b_identity"] == "not-passed"
    assert "annotation_generation_independence" not in accepted["blocking_gates"]


def test_candidate_acceptance_rejects_overclaimed_license(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    qualification = json.loads((run_dir / "source_qualification.json").read_text(encoding="utf-8"))
    qualification["gates"]["raw_file_license"] = "passed"
    (run_dir / "source_qualification.json").write_text(json.dumps(qualification), encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_source_tampered_after_independent_check(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    (tmp_path / "archive.sra").write_bytes(b"tampered-after-independent")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_synchronized_gate_overclaim(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    qualification_path = run_dir / "source_qualification.json"
    qualification = json.loads(qualification_path.read_text(encoding="utf-8"))
    qualification["gates"]["unmodified_strain-b_identity"] = "passed"
    qualification_path.write_text(json.dumps(qualification, sort_keys=True), encoding="utf-8")
    shutil.copy2(qualification_path, peer / "source_qualification.json")
    run_manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    run_manifest["artifacts"]["source_qualification.json"] = {"size_bytes": qualification_path.stat().st_size, "sha256": hashlib.sha256(qualification_path.read_bytes()).hexdigest()}
    (run_dir / "run_manifest.json").write_text(json.dumps(run_manifest, sort_keys=True), encoding="utf-8")
    shutil.copy2(run_dir / "run_manifest.json", peer / "run_manifest.json")
    independent_value = json.loads(independent.read_text(encoding="utf-8"))
    independent_value["recomputed_gates"] = qualification["gates"]
    independent_value["verified_artifacts"]["source_qualification.json"] = run_manifest["artifacts"]["source_qualification.json"]["sha256"]
    independent.write_text(json.dumps(independent_value, sort_keys=True), encoding="utf-8")
    shutil.copy2(independent, peer / "verification/independent_check.json")
    repeat_value = json.loads(repeat.read_text(encoding="utf-8"))
    for relative in REPEATABLE_FILES:
        repeat_value["file_sha256"][relative] = hashlib.sha256((run_dir / relative).read_bytes()).hexdigest()
    repeat.write_text(json.dumps(repeat_value, sort_keys=True), encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_tampered_protected_artifact(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    (tmp_path / "slice0/evidence.json").write_text("tampered", encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_receipt_input_hash_mismatch(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    receipt_path = peer / "verification/execution_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["manifest_sha256"] = "0" * 64
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    repeat_value = json.loads(repeat.read_text(encoding="utf-8"))
    repeat_value["execution_receipts"][peer.name]["identity"] = _identity(receipt_path)
    repeat.write_text(json.dumps(repeat_value), encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_extra_run_artifact(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    (run_dir / "candidate.tsv").write_text("out-of-scope", encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_forged_workflow_argv(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    receipt_path = peer / "verification/execution_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["workflow_steps"][0]["command"] = receipt["workflow_steps"][0]["command"][:2]
    receipt["command"] = receipt["workflow_steps"][0]["command"]
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    repeat_value = json.loads(repeat.read_text(encoding="utf-8"))
    repeat_value["execution_receipts"][peer.name]["identity"] = _identity(receipt_path)
    repeat.write_text(json.dumps(repeat_value), encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)


def test_candidate_acceptance_rejects_out_of_order_workflow_times(tmp_path: Path) -> None:
    manifest, reference, root = _fixture(tmp_path)
    run_dir = tmp_path / "run"
    qualify_independent_transcript_source(manifest, reference, run_dir, root)
    peer = tmp_path / "peer"
    independent, tests, repeat = _verification(run_dir, peer, tmp_path, manifest, reference, root)
    receipt_path = peer / "verification/execution_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["workflow_steps"][1]["started_at_utc"] = "2026-07-14T23:59:59+00:00"
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    repeat_value = json.loads(repeat.read_text(encoding="utf-8"))
    repeat_value["execution_receipts"][peer.name]["identity"] = _identity(receipt_path)
    repeat.write_text(json.dumps(repeat_value), encoding="utf-8")
    with pytest.raises(ContractError):
        create_independent_source_acceptance(run_dir, peer, independent, tests, repeat, root, manifest, reference)
