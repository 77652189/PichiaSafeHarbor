from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.slice0a import create_slice0a_acceptance


SOURCE_IDS = [
    "strain-b_insdc_gca001746955_1",
    "strain-b_refseq_gcf000027005_1",
    "strain-c_genbank_gca000223565_1",
    "ensembl_fungi_release63_strain-b",
    "kegg_ppa_metadata",
    "pichiagenome_2016_curation",
    "biocyc_picpa",
]
SEQUENCES = ["NC_012963.1", "NC_012964.1", "NC_012965.1", "NC_012966.1"]
STATISTICS = {
    "gene_count": 1,
    "transcript_count": 1,
    "orf_cds_count": 1,
    "ncrna_count": 0,
    "trna_count": 0,
    "rrna_count": 0,
    "partial_gene_count": 0,
    "start_range_gene_count": 0,
    "end_range_gene_count": 0,
    "sequence_count": 1,
}


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def _identity(path: Path) -> dict[str, int | str]:
    return {"size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _source(source_id: str) -> dict:
    unavailable = source_id in {"pichiagenome_2016_curation", "biocyc_picpa"}
    independence_group = "strain-b_gca000027005_insdc" if source_id in {
        "strain-b_refseq_gcf000027005_1",
        "ensembl_fungi_release63_strain-b",
        "kegg_ppa_metadata",
    } else source_id
    relationship = {
        "ensembl_fungi_release63_strain-b": "mirror_reformat_of_old_insdc",
        "kegg_ppa_metadata": "functional_metadata_mirror_of_old_annotation",
    }.get(source_id, "candidate")
    return {
        "source_id": source_id,
        "source_version": "fixture-v1",
        "source_url": "https://example.invalid/source",
        "availability": "unavailable" if unavailable else "available",
        "coordinate_space": "GCA_001746955.1",
        "license": {"status": "fixture-license", "url": "https://example.invalid/license"},
        "files": {} if unavailable else {"annotation": {"path": "fixture.gff", "size_bytes": 1, "sha256": "a" * 64}},
        "upstream": {"independence_group": independence_group, "relationship": relationship},
        "experimental_support": ["RNA-seq", "proteomics"] if source_id == "pichiagenome_2016_curation" else [],
    }


def _mapping() -> dict:
    records = []
    statuses = ["conflict", "unmappable", "uncertain"]
    for sequence_index, seqid in enumerate(SEQUENCES):
        for position_index, position in enumerate(("first", "middle", "last")):
            window_status = "uncertain" if sequence_index == 0 and position_index == 0 else "preserved"
            for gene_index in range(3):
                status = statuses.pop(0) if statuses else "consistent"
                records.append(
                    {
                        "source_seqid": seqid,
                        "window_id": f"{seqid}:{position}",
                        "window_position": position,
                        "mapping_status": status,
                        "order_status": window_status,
                        "adjacency_status": "unavailable",
                    }
                )
    mapping_counts = dict(sorted(Counter(item["mapping_status"] for item in records).items()))
    return {
        "records": records,
        "summary": {
            "mapping_status_counts": mapping_counts,
            "source_nuclear_sequences_covered": SEQUENCES,
            "order_status_counts": {"preserved": 11, "uncertain": 1},
            "adjacency_status_counts": {"unavailable": 12},
        },
    }


def _fixture_run(tmp_path: Path) -> tuple[Path, Path, Path]:
    run_dir = tmp_path / "run"
    sources = [_source(source_id) for source_id in SOURCE_IDS]
    _write_json(
        run_dir / "annotation_sources.json",
        {
            "primary_coordinate_space": "GCA_001746955.1",
            "sources": sources,
            "acquisition_evidence": {
                "probes": [
                    {
                        "probe_id": "pichiagenome_http_redirect",
                        "target": "http://www.pichiagenome.org",
                        "method": "curl -L --max-time 30",
                        "result": "HTTP 502",
                        "availability": "unavailable",
                    }
                ]
            },
        },
    )
    qualifications = []
    for source_id in SOURCE_IDS:
        qualifications.append(
            {
                "source_id": source_id,
                "qualification_status": "unavailable" if source_id in {"pichiagenome_2016_curation", "biocyc_picpa"} else "unqualified",
                "statistics": None if source_id in {"pichiagenome_2016_curation", "biocyc_picpa"} else STATISTICS,
            }
        )
    _write_json(
        run_dir / "annotation_qualification.json",
        {
            "sources": qualifications,
            "mirror_evidence": {
                "ensembl_vs_old_refseq": {"shared_gene_id_count": 5040, "exact_boundary_match_count": 5040},
                "kegg_vs_old_refseq": {"shared_gene_id_count": 5040, "exact_coordinate_match_count": 5040},
            },
        },
    )
    mapping = _mapping()
    _write_json(run_dir / "annotation_mapping_probe.json", mapping)
    (run_dir / "annotation_qualification.tsv").write_text("fixture\n", encoding="utf-8")
    (run_dir / "annotation_mapping_probe.tsv").write_text("fixture\n", encoding="utf-8")
    (run_dir / "annotation_qualification_report.md").write_text(
        "The 28 high-confidence intervals are not accepted.\nNo thresholds, candidate windows, or formal candidates are generated.\n",
        encoding="utf-8",
    )
    (run_dir / "source_selection_recommendation.md").write_text(
        "Do not start it without a new ADR.\nScientific status: blocked.\n".replace("Do not", "do not"),
        encoding="utf-8",
    )
    artifact_names = [
        "annotation_sources.json",
        "annotation_qualification.json",
        "annotation_qualification.tsv",
        "annotation_mapping_probe.json",
        "annotation_mapping_probe.tsv",
        "annotation_qualification_report.md",
        "source_selection_recommendation.md",
    ]
    artifacts = {name: _identity(run_dir / name) for name in artifact_names}
    _write_json(
        run_dir / "run_manifest.json",
        {
            "run_id": "slice0a-0123456789abcdef",
            "implementation_sha256": "a" * 64,
            "primary_coordinate_space": "GCA_001746955.1",
            "scientific_acceptance_status": "blocked",
            "artifacts": artifacts,
        },
    )
    independent = run_dir / "verification/independent_check.json"
    _write_json(
        independent,
        {
            "schema_version": 1,
            "verification_type": "slice0a_independent_qualification",
            "status": "passed",
            "run_id": "slice0a-0123456789abcdef",
            "source_count": 7,
            "source_statistics": {source_id: STATISTICS for source_id in SOURCE_IDS[:5]},
            "probe_record_count": 36,
            "mapping_status_counts": mapping["summary"]["mapping_status_counts"],
            "order_status_counts": mapping["summary"]["order_status_counts"],
            "adjacency_status_counts": mapping["summary"]["adjacency_status_counts"],
            "source_nuclear_sequences_covered": SEQUENCES,
            "mirror_checks": {"ensembl_shared": 5040, "ensembl_exact": 5040, "kegg_shared": 5040, "kegg_exact": 5040},
            "verified_artifacts": {name: identity["sha256"] for name, identity in artifacts.items()},
        },
    )
    tests = run_dir / "verification/test_evidence.json"
    _write_json(
        tests,
        {
            "schema_version": 1,
            "evidence_type": "automated_tests",
            "status": "passed",
            "command": ["python", "-m", "pytest", "-q"],
            "passed_count": 1,
            "implementation_sha256": "a" * 64,
        },
    )
    return run_dir, independent, tests


def _refresh_artifact(run_dir: Path, independent: Path, name: str) -> None:
    run = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    run["artifacts"][name] = _identity(run_dir / name)
    _write_json(run_dir / "run_manifest.json", run)
    evidence = json.loads(independent.read_text(encoding="utf-8"))
    evidence["verified_artifacts"][name] = run["artifacts"][name]["sha256"]
    _write_json(independent, evidence)


def test_slice0a_acceptance_binds_independent_source_statistics(tmp_path: Path) -> None:
    run_dir, independent, tests = _fixture_run(tmp_path)
    assert create_slice0a_acceptance(run_dir, independent, tests)["verification_status"] == "passed"
    evidence = json.loads(independent.read_text(encoding="utf-8"))
    evidence["source_statistics"][SOURCE_IDS[0]]["gene_count"] = 2
    _write_json(independent, evidence)
    with pytest.raises(ContractError, match="source statistics mismatch"):
        create_slice0a_acceptance(run_dir, independent, tests)


@pytest.mark.parametrize("failure", ["missing_source", "missing_license", "missing_mapping_class", "missing_window"])
def test_slice0a_acceptance_rejects_incomplete_completion_evidence(tmp_path: Path, failure: str) -> None:
    run_dir, independent, tests = _fixture_run(tmp_path)
    if failure in {"missing_source", "missing_license"}:
        path = run_dir / "annotation_sources.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        if failure == "missing_source":
            value["sources"] = value["sources"][:-1]
        else:
            value["sources"][0]["license"] = {}
        _write_json(path, value)
        _refresh_artifact(run_dir, independent, path.name)
    else:
        path = run_dir / "annotation_mapping_probe.json"
        value = json.loads(path.read_text(encoding="utf-8"))
        if failure == "missing_mapping_class":
            for record in value["records"]:
                if record["mapping_status"] == "conflict":
                    record["mapping_status"] = "consistent"
            value["summary"]["mapping_status_counts"] = dict(sorted(Counter(item["mapping_status"] for item in value["records"]).items()))
        else:
            for record in value["records"]:
                if record["source_seqid"] == SEQUENCES[0] and record["window_position"] == "last":
                    record["window_position"] = "middle"
        _write_json(path, value)
        _refresh_artifact(run_dir, independent, path.name)
        evidence = json.loads(independent.read_text(encoding="utf-8"))
        evidence["mapping_status_counts"] = value["summary"]["mapping_status_counts"]
        _write_json(independent, evidence)
    with pytest.raises(ContractError):
        create_slice0a_acceptance(run_dir, independent, tests)
