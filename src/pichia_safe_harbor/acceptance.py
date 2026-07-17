from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ContractError(f"required acceptance evidence is missing: {path}")
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ContractError(f"acceptance evidence must be a JSON object: {path}")
    return value


def create_acceptance_manifest(
    run_dir: Path,
    independent_check_path: Path,
    test_evidence_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    run_manifest_path = run_dir / "run_manifest.json"
    run_manifest = _read_json(run_manifest_path)
    if run_manifest.get("schema_version") != 2:
        raise ContractError("acceptance requires run manifest schema_version 2")
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

    statistics = _read_json(run_dir / "statistics.json")
    independent = _read_json(independent_check_path)
    if independent.get("schema_version") != 1 or independent.get(
        "verification_type"
    ) != "independent_interval_recalculation":
        raise ContractError("unsupported independent verification evidence schema")
    if independent.get("status") != "passed":
        raise ContractError("independent verification did not pass")
    if independent.get("run_id") != run_manifest["run_id"]:
        raise ContractError("independent verification run_id mismatch")
    independent_artifacts = independent.get("verified_artifacts", {})
    for name in ("intergenic_regions.tsv", "terminal_regions.tsv"):
        if independent_artifacts.get(name) != run_manifest["artifacts"][name]["sha256"]:
            raise ContractError(
                f"independent verification artifact identity mismatch: {name}"
            )
    if (
        independent.get("verified_intergenic_region_count")
        != statistics["intergenic_region_count"]
    ):
        raise ContractError("independent verification intergenic count mismatch")
    expected_nuclear = sorted(
        seqid
        for seqid, values in statistics["sequence_statistics"].items()
        if not values["excluded_from_nuclear_statistics"]
    )
    if independent.get("nuclear_sequences") != expected_nuclear:
        raise ContractError("independent verification nuclear sequence set mismatch")
    coverage_checks = independent.get("coverage_checks", {})
    if sorted(coverage_checks) != expected_nuclear or not all(
        values.get("conserved") is True for values in coverage_checks.values()
    ):
        raise ContractError("independent verification coverage conservation mismatch")

    test_evidence = _read_json(test_evidence_path)
    if test_evidence.get("schema_version") != 1 or test_evidence.get(
        "evidence_type"
    ) != "automated_tests":
        raise ContractError("unsupported automated test evidence schema")
    if test_evidence.get("status") != "passed":
        raise ContractError("automated tests did not pass")
    if test_evidence.get("command") != ["python", "-m", "pytest", "-q"]:
        raise ContractError("automated test evidence command mismatch")
    if not isinstance(test_evidence.get("passed_count"), int) or test_evidence[
        "passed_count"
    ] <= 0:
        raise ContractError("automated test evidence passed_count must be positive")
    if (
        test_evidence.get("implementation_sha256")
        != run_manifest["software"]["implementation_sha256"]
    ):
        raise ContractError("test evidence implementation hash mismatch")

    endpoint_summary = {
        seqid: values["endpoint_completeness"]
        for seqid, values in sorted(statistics["sequence_statistics"].items())
        if not values["excluded_from_nuclear_statistics"]
    }
    blockers = [
        "systematic_annotation_boundary_uncertainty_requires_authoritative_review",
        "authoritative_scientific_acceptance_pending",
    ]
    result = {
        "schema_version": 1,
        "acceptance_version": "slice0-adr0003-v1",
        "run_id": run_manifest["run_id"],
        "execution_status": "complete",
        "verification_status": "passed",
        "scientific_acceptance_status": "blocked",
        "scientific_acceptance_blockers": blockers,
        "inputs": run_manifest["inputs"],
        "reference_identity": {
            "target_strain": run_manifest["target_strain"],
            "primary_reference_strain": run_manifest["primary_reference_strain"],
            "primary_assembly": run_manifest["primary_assembly"],
            "secondary_reference_strain": run_manifest["secondary_reference_strain"],
            "secondary_assembly": run_manifest["secondary_assembly"],
            "exact_target_strain_coordinates": run_manifest["exact_target_strain_coordinates"],
        },
        "sequence_name_mapping": run_manifest["sequence_name_mapping"],
        "sequence_classes": run_manifest["sequence_classes"],
        "implementation": run_manifest["software"],
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
            "independent_interval_check": {
                "size_bytes": independent_check_path.stat().st_size,
                "sha256": sha256_file(independent_check_path),
                "verified_intergenic_region_count": independent[
                    "verified_intergenic_region_count"
                ],
            },
        },
        "completeness_summary": {
            "nuclear_sequence_endpoints": endpoint_summary,
            "annotation_boundaries": statistics["annotation_boundary_summary"],
            "all_intergenic_region_count": statistics["intergenic_region_count"],
            "high_confidence_intergenic": statistics["high_confidence_intergenic"],
        },
    }
    destination = output_path or (run_dir / "acceptance_manifest.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result
