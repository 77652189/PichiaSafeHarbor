from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .candidates import RULE_VERSION, classify_candidate_windows
from .errors import ContractError
from .io_utils import sha256_file, write_json, write_tsv
from .models import ExcludedRegion, IntergenicRegion
from .parsers import read_fasta_sequences
from .pipeline import _implementation_hash


def _write_fasta(path: Path, records: Iterable[dict[str, Any]]) -> None:
    lines: list[str] = []
    for record in records:
        lines.append(
            f">{record['candidate_id']} seqid={record['seqid']} "
            f"start={record['start']} end={record['end']} length={record['length']}"
        )
        sequence = record.get("sequence", "")
        for offset in range(0, len(sequence), 70):
            lines.append(sequence[offset : offset + 70])
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _load_intergenic_regions(baseline_run_dir: Path) -> list[IntergenicRegion]:
    path = baseline_run_dir / "intergenic_regions.json"
    if not path.is_file():
        raise ContractError(f"baseline run is missing intergenic_regions.json: {path}")
    records = json.loads(path.read_text(encoding="utf-8"))
    return [
        IntergenicRegion(
            region_id=record["region_id"],
            seqid=record["seqid"],
            start=record["start"],
            end=record["end"],
            length=record["length"],
            left_entity_ids=tuple(record["left_entity_ids"]),
            right_entity_ids=tuple(record["right_entity_ids"]),
            left_strand=record["left_strand"],
            right_strand=record["right_strand"],
            orientation=record["orientation"],
            gap_length_bp=record["gap_length_bp"],
            left_boundary_confidence=record["left_boundary_confidence"],
            right_boundary_confidence=record["right_boundary_confidence"],
            coordinate_system=record.get("coordinate_system", "0-based-half-open"),
        )
        for record in records
    ]


def _run_id(parent_run_id: str, rule_params: dict[str, Any], implementation_hash: str) -> str:
    material = json.dumps(
        {
            "parent_run_id": parent_run_id,
            "rule_params": rule_params,
            "rule_version": RULE_VERSION,
            "implementation_sha256": implementation_hash,
        },
        sort_keys=True,
    ).encode("utf-8")
    return "slice1-" + hashlib.sha256(material).hexdigest()[:16]


def run_slice1(
    baseline_run_dir: Path,
    reference: dict[str, Any],
    fasta_path: Path,
    buffer_distance_bp: int,
    min_candidate_window_bp: int,
    output_dir: Path,
    long_interval_percentile: float = 0.95,
) -> dict[str, Any]:
    if output_dir.exists():
        raise ContractError(f"output directory already exists: {output_dir}")
    manifest_path = baseline_run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise ContractError(f"baseline run is missing run_manifest.json: {manifest_path}")
    baseline_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if baseline_manifest.get("execution_status") != "complete":
        raise ContractError(
            "baseline run execution_status is not complete: "
            f"{baseline_manifest.get('execution_status')!r}"
        )

    baseline_fasta_sha256 = baseline_manifest.get("inputs", {}).get("fasta", {}).get("sha256")
    fasta_sha256 = sha256_file(fasta_path)
    if baseline_fasta_sha256 is not None and fasta_sha256 != baseline_fasta_sha256:
        raise ContractError(
            "fasta file does not match the one used for the parent baseline run: "
            f"expected sha256={baseline_fasta_sha256}, got {fasta_sha256} ({fasta_path})"
        )
    sequences = read_fasta_sequences(fasta_path)

    intergenic_regions = _load_intergenic_regions(baseline_run_dir)

    rule_params = {
        "buffer_distance_bp": buffer_distance_bp,
        "min_candidate_window_bp": min_candidate_window_bp,
        "long_interval_percentile": long_interval_percentile,
    }
    windows, excluded, summary = classify_candidate_windows(
        intergenic_regions, buffer_distance_bp, min_candidate_window_bp, long_interval_percentile
    )
    windows = [
        replace(window, sequence=sequences[window.seqid][window.start : window.end])
        for window in windows
    ]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        window_records = [item.to_dict() for item in windows]
        excluded_records = [item.to_dict() for item in excluded]
        write_json(temp_dir / "candidate_windows.json", window_records)
        write_tsv(
            temp_dir / "candidate_windows.tsv",
            window_records,
            [
                "candidate_id",
                "parent_region_id",
                "seqid",
                "start",
                "end",
                "length",
                "sequence",
                "left_entity_ids",
                "right_entity_ids",
                "orientation",
                "left_boundary_confidence",
                "right_boundary_confidence",
                "evidence_level",
                "structural_tier",
                "rule_flags",
                "coordinate_system",
            ],
        )
        _write_fasta(temp_dir / "candidate_windows.fasta", window_records)
        write_json(temp_dir / "excluded_regions.json", excluded_records)
        write_tsv(
            temp_dir / "excluded_regions.tsv",
            excluded_records,
            list(ExcludedRegion.__dataclass_fields__),
        )
        write_json(temp_dir / "statistics.json", summary)
        report = _render_report(reference, rule_params, summary)
        (temp_dir / "candidate_windows_report.md").write_text(report, encoding="utf-8")

        artifact_names = [
            "candidate_windows.json",
            "candidate_windows.tsv",
            "candidate_windows.fasta",
            "excluded_regions.json",
            "excluded_regions.tsv",
            "statistics.json",
            "candidate_windows_report.md",
        ]
        implementation_hash = _implementation_hash()
        run_manifest = {
            "schema_version": 1,
            "run_id": _run_id(baseline_manifest["run_id"], rule_params, implementation_hash),
            "parent_run_id": baseline_manifest["run_id"],
            "execution_status": "complete",
            "verification_status": "not_run",
            "scientific_acceptance_status": "blocked",
            "scientific_acceptance_blockers": [
                "rule_parameters_are_illustrative_placeholders_pending_threshold_adr",
                "input_annotation_is_strain-b_genbank_2016_release_not_an_adr_qualified_boundary_track",
                "no_target_strain_native_sequencing_data_used",
                "collinearity_with_strain-c_not_yet_computed",
                "acceptance_manifest_not_yet_generated",
            ],
            "run_purpose": "engine_readiness_test_not_scientific_output",
            "software": {
                "name": "pichia-safe-harbor",
                "version": __version__,
                "implementation_sha256": implementation_hash,
            },
            "rule_version": RULE_VERSION,
            "rule_params": rule_params,
            "target_strain": "Strain-T",
            "primary_reference_strain": "Strain-B",
            "primary_assembly": "GCA_001746955.1",
            "secondary_reference_strain": "Strain-C",
            "secondary_assembly": "GCA_000223565.1",
            "strain_applicability": "close-strain proxy",
            "exact_target_strain_coordinates": False,
            "analysis_scope": "Slice 1 candidate windows over a Slice 0 baseline; placeholder rule parameters",
            "coordinate_system": "0-based-half-open",
            "known_limitations": reference.get("known_limitations", []),
            "missing_risk_tracks": reference.get("missing_risk_tracks", []),
            "candidate_window_count": summary["candidate_window_count"],
            "excluded_region_count": summary["excluded_region_count"],
            "artifacts": {
                name: {
                    "sha256": sha256_file(temp_dir / name),
                    "size_bytes": (temp_dir / name).stat().st_size,
                }
                for name in artifact_names
            },
        }
        write_json(temp_dir / "run_manifest.json", run_manifest)
        os.replace(temp_dir, output_dir)
        return run_manifest
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


def create_slice1_acceptance(
    run_dir: Path,
    independent_check_path: Path,
    test_evidence_path: Path,
) -> dict[str, Any]:
    run_manifest_path = run_dir / "run_manifest.json"
    if not run_manifest_path.is_file():
        raise ContractError(f"Slice 1 run manifest is missing: {run_manifest_path}")
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    if run_manifest.get("schema_version") != 1:
        raise ContractError("acceptance requires Slice 1 run manifest schema_version 1")
    if run_manifest.get("execution_status") != "complete":
        raise ContractError("acceptance requires execution_status=complete")

    verified_artifacts: dict[str, dict[str, Any]] = {}
    for name, expected in sorted(run_manifest.get("artifacts", {}).items()):
        path = run_dir / name
        if not path.is_file():
            raise ContractError(f"run artifact is missing: {path}")
        actual = {"size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
        if actual != expected:
            raise ContractError(f"run artifact identity mismatch: {name}")
        verified_artifacts[name] = actual

    statistics = json.loads((run_dir / "statistics.json").read_text(encoding="utf-8"))

    independent = json.loads(independent_check_path.read_text(encoding="utf-8"))
    if (
        independent.get("schema_version") != 1
        or independent.get("verification_type") != "slice1_independent_recalculation"
    ):
        raise ContractError("unsupported independent verification evidence schema")
    if independent.get("status") != "passed":
        raise ContractError("independent verification did not pass")
    if independent.get("run_id") != run_manifest["run_id"]:
        raise ContractError("independent verification run_id mismatch")
    if independent.get("verified_artifacts") != {
        name: identity["sha256"] for name, identity in run_manifest["artifacts"].items()
    }:
        raise ContractError("independent verification artifact identity mismatch")
    if independent.get("verified_candidate_window_count") != statistics["candidate_window_count"]:
        raise ContractError("independent verification candidate window count mismatch")
    if independent.get("verified_excluded_region_count") != statistics["excluded_region_count"]:
        raise ContractError("independent verification excluded region count mismatch")
    if independent.get("verified_excluded_by_rule") != statistics["excluded_by_rule"]:
        raise ContractError("independent verification excluded-by-rule mismatch")

    test_evidence = json.loads(test_evidence_path.read_text(encoding="utf-8"))
    if (
        test_evidence.get("schema_version") != 1
        or test_evidence.get("evidence_type") != "automated_tests"
    ):
        raise ContractError("unsupported automated test evidence schema")
    if test_evidence.get("status") != "passed":
        raise ContractError("automated tests did not pass")
    if test_evidence.get("command") != ["python", "-m", "pytest", "-q"]:
        raise ContractError("automated test evidence command mismatch")
    if not isinstance(test_evidence.get("passed_count"), int) or test_evidence["passed_count"] <= 0:
        raise ContractError("automated test evidence passed_count must be positive")
    if (
        test_evidence.get("implementation_sha256")
        != run_manifest["software"]["implementation_sha256"]
    ):
        raise ContractError("test evidence implementation hash mismatch")

    result = {
        "schema_version": 1,
        "acceptance_version": "slice1-adr0011-v1",
        "run_id": run_manifest["run_id"],
        "parent_run_id": run_manifest["parent_run_id"],
        "execution_status": "complete",
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
        "scientific_acceptance_blockers": run_manifest["scientific_acceptance_blockers"],
        "run_purpose": run_manifest["run_purpose"],
        "rule_version": run_manifest["rule_version"],
        "rule_params": run_manifest["rule_params"],
        "run_manifest": {
            "size_bytes": run_manifest_path.stat().st_size,
            "sha256": sha256_file(run_manifest_path),
        },
        "artifacts": verified_artifacts,
        "verification_evidence": {
            "automated_tests": {
                "size_bytes": test_evidence_path.stat().st_size,
                "sha256": sha256_file(test_evidence_path),
                "command": test_evidence["command"],
                "passed_count": test_evidence["passed_count"],
            },
            "independent_recalculation": {
                "size_bytes": independent_check_path.stat().st_size,
                "sha256": sha256_file(independent_check_path),
                "verified_candidate_window_count": independent["verified_candidate_window_count"],
                "verified_excluded_region_count": independent["verified_excluded_region_count"],
            },
        },
        "candidate_summary": {
            "input_intergenic_region_count": statistics["input_intergenic_region_count"],
            "candidate_window_count": statistics["candidate_window_count"],
            "excluded_region_count": statistics["excluded_region_count"],
            "excluded_by_rule": statistics["excluded_by_rule"],
            "candidate_structural_tier_counts": statistics["candidate_structural_tier_counts"],
        },
    }
    destination = run_dir / "acceptance_manifest.json"
    destination.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


def _render_report(reference: dict[str, Any], rule_params: dict[str, Any], summary: dict[str, Any]) -> str:
    lines = [
        "# Strain-B Slice 1 candidate safe-harbor windows (ENGINE READINESS TEST -- NOT SCIENTIFIC OUTPUT)",
        "",
        "## Status",
        "",
        "- run_purpose: engine_readiness_test_not_scientific_output",
        "- scientific_acceptance_status: blocked",
        "- Target strain: Strain-T (this run is a close-strain Strain-B proxy; exact_target_strain_coordinates=false)",
        "- Input annotation: Strain-B GenBank 2016-09-21 release; not an ADR-qualified boundary track",
        "- Rule parameters below are ILLUSTRATIVE PLACEHOLDERS, not an ADR-frozen threshold",
        "",
        "## Rule parameters (placeholder, unvalidated)",
        "",
        f"- rule_version: {RULE_VERSION}",
        f"- buffer_distance_bp: {rule_params['buffer_distance_bp']}",
        f"- min_candidate_window_bp: {rule_params['min_candidate_window_bp']}",
        f"- long_interval_percentile: {rule_params['long_interval_percentile']}",
        "",
        "## Results",
        "",
        f"- Input intergenic regions: {summary['input_intergenic_region_count']}",
        f"- Candidate windows produced: {summary['candidate_window_count']}",
        f"- Excluded regions: {summary['excluded_region_count']}",
        f"- Excluded by rule: {json.dumps(summary['excluded_by_rule'], sort_keys=True)}",
        f"- Candidate structural tiers: {json.dumps(summary['candidate_structural_tier_counts'], sort_keys=True)}",
        f"- Long-interval length threshold (bp): {summary['long_interval_length_threshold_bp']}",
        "",
        "## Unresolved global evidence gaps (apply to every candidate window in this run)",
        "",
    ]
    for item in reference.get("missing_risk_tracks", []):
        lines.append(f"- {item}: unavailable")
    lines.extend(
        [
            "- collinearity_status: unavailable (Slice 2 Strain-C concordance not yet implemented)",
            "- target_strain_specific_variation: unavailable (no Strain-T native sequencing data used)",
            "",
            "## Known reference limitations",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in reference.get("known_limitations", []))
    lines.append("")
    return "\n".join(lines) + "\n"
