from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable

from . import __version__
from .analysis import (
    build_inventory,
    map_annotation_sequences,
    normalize_entities,
    validate_coordinates,
    validate_identity,
)
from .errors import ContractError
from .io_utils import sha256_file, write_json
from .models import FunctionalEntity, IntergenicRegion, TerminalRegion
from .parsers import feature_type_counts, read_fasta_index, read_gff3
from .reference import verify_reference_paths


def _write_tsv(path: Path, records: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    key: ",".join(value) if isinstance(value, (list, tuple)) else value
                    for key, value in record.items()
                }
            )


def _implementation_hash() -> str:
    digest = hashlib.sha256()
    package_dir = Path(__file__).parent
    for path in sorted(package_dir.glob("*.py"), key=lambda item: item.name):
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _run_id(
    reference: dict[str, Any],
    implementation_hash: str,
    inputs: dict[str, dict[str, Any]],
    applied_mapping: dict[str, str],
) -> str:
    material = json.dumps(
        {
            "reference_identity": {
                "key": reference["key"],
                "assembly_accession": reference["assembly_accession"],
                "annotation_accession": reference["annotation_accession"],
                "manifest_version": reference["manifest_version"],
            },
            "actual_inputs": {
                role: {
                    "size_bytes": values["size_bytes"],
                    "sha256": values["sha256"],
                }
                for role, values in inputs.items()
            },
            "sequence_name_mapping": applied_mapping,
            "sequence_classes": reference["sequence_classes"],
            "sequence_endpoints": reference["sequence_endpoints"],
            "software_version": __version__,
            "implementation_sha256": implementation_hash,
            "coordinate_system": "0-based-half-open",
            "boundary_confidence_policy": "adr-0003-v1",
        },
        sort_keys=True,
    ).encode("utf-8")
    return "slice0-" + hashlib.sha256(material).hexdigest()[:16]


def run_baseline(reference: dict[str, Any], fasta_path: Path, annotation_path: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists():
        raise ContractError(f"output directory already exists: {output_dir}")
    verified_paths = verify_reference_paths(reference, fasta_path, annotation_path)
    inputs = {
        role: {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for role, path in verified_paths.items()
    }
    fasta = read_fasta_index(fasta_path)
    features, directives, sequence_regions = read_gff3(annotation_path)
    validate_identity(reference, directives)
    mapped_features, applied_mapping = map_annotation_sequences(
        features, fasta, reference.get("sequence_name_map", {})
    )
    validate_coordinates(mapped_features, fasta)
    for seqid, (_start, end) in sequence_regions.items():
        mapped_seqid = applied_mapping.get(seqid, seqid)
        if mapped_seqid in fasta and end != fasta[mapped_seqid].length:
            raise ContractError(
                f"GFF sequence-region length mismatch for {mapped_seqid}: {end} != {fasta[mapped_seqid].length}"
            )
    entities, diagnostics = normalize_entities(mapped_features)
    intergenic, terminal, statistics = build_inventory(
        fasta,
        entities,
        reference["sequence_classes"],
        reference["sequence_endpoints"],
    )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".", dir=output_dir.parent))
    try:
        intergenic_records = [record.to_dict() for record in intergenic]
        terminal_records = [record.to_dict() for record in terminal]
        entity_records = [record.to_dict() for record in entities]
        write_json(temp_dir / "functional_entities.json", entity_records)
        _write_tsv(
            temp_dir / "functional_entities.tsv",
            entity_records,
            list(FunctionalEntity.__dataclass_fields__),
        )
        write_json(temp_dir / "intergenic_regions.json", intergenic_records)
        _write_tsv(
            temp_dir / "intergenic_regions.tsv",
            intergenic_records,
            list(IntergenicRegion.__dataclass_fields__),
        )
        write_json(temp_dir / "terminal_regions.json", terminal_records)
        _write_tsv(
            temp_dir / "terminal_regions.tsv",
            terminal_records,
            list(TerminalRegion.__dataclass_fields__),
        )
        statistics["feature_type_counts"] = feature_type_counts(mapped_features)
        statistics["diagnostic_counts"] = {
            key: len(value) if isinstance(value, list) else sum(value.values())
            for key, value in diagnostics.to_dict().items()
        }
        write_json(temp_dir / "statistics.json", statistics)
        write_json(temp_dir / "diagnostics.json", diagnostics.to_dict())
        limitations = reference["known_limitations"]
        report = _render_report(reference, statistics, limitations)
        (temp_dir / "baseline_report.md").write_text(report, encoding="utf-8")

        artifact_names = [
            "functional_entities.json",
            "functional_entities.tsv",
            "intergenic_regions.json",
            "intergenic_regions.tsv",
            "terminal_regions.json",
            "terminal_regions.tsv",
            "statistics.json",
            "diagnostics.json",
            "baseline_report.md",
        ]
        implementation_hash = _implementation_hash()
        run_manifest = {
            "schema_version": 2,
            "run_id": _run_id(reference, implementation_hash, inputs, applied_mapping),
            "execution_status": "complete",
            "verification_status": "not_run",
            "scientific_acceptance_status": "blocked",
            "scientific_acceptance_blockers": [
                "systematic_annotation_boundary_uncertainty_requires_authoritative_review",
                "acceptance_manifest_not_yet_generated",
            ],
            "software": {
                "name": "pichia-safe-harbor",
                "version": __version__,
                "implementation_sha256": implementation_hash,
            },
            "target_strain": "Strain-T",
            "primary_reference_strain": "Strain-B",
            "primary_assembly": "GCA_001746955.1",
            "secondary_reference_strain": "Strain-C",
            "secondary_assembly": "GCA_000223565.1",
            "strain_applicability": "close-strain proxy",
            "exact_target_strain_coordinates": False,
            "analysis_scope": "Slice 0 raw intergenic baseline; no candidate windows or thresholds",
            "coordinate_system": "0-based-half-open",
            "reference_manifest_version": reference["manifest_version"],
            "annotation_accession": reference["annotation_accession"],
            "inputs": inputs,
            "sequence_name_mapping": applied_mapping,
            "sequence_classes": reference["sequence_classes"],
            "sequence_endpoints": reference["sequence_endpoints"],
            "known_limitations": limitations,
            "missing_risk_tracks": reference["missing_risk_tracks"],
            "artifacts": {
                name: {"sha256": sha256_file(temp_dir / name), "size_bytes": (temp_dir / name).stat().st_size}
                for name in artifact_names
            },
        }
        write_json(temp_dir / "run_manifest.json", run_manifest)
        os.replace(temp_dir, output_dir)
        return run_manifest
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


def _render_report(reference: dict[str, Any], statistics: dict[str, Any], limitations: list[str]) -> str:
    lines = [
        "# Strain-B Slice 0 complete-genome intergenic baseline",
        "",
        "## Applicability",
        "",
        "- Target strain: Strain-T",
        "- Primary coordinate reference: Strain-B GCA_001746955.1",
        "- Strain applicability: close-strain proxy",
        "- Exact Strain-T coordinates: false",
        "- Scope: raw intergenic facts only; no candidate windows or frozen thresholds",
        "",
        "## Baseline summary",
        "",
        f"- Nuclear chromosomes analyzed: {statistics['nuclear_chromosome_count']}",
        f"- Raw intergenic regions: {statistics['intergenic_region_count']}",
        f"- Orientation counts: {json.dumps(statistics['orientation_counts'], sort_keys=True)}",
        f"- High-confidence boundary subset: {statistics['high_confidence_intergenic']['count']}",
        f"- Partial functional entities: {statistics['annotation_boundary_summary']['partial_entity_count']} / {statistics['annotation_boundary_summary']['functional_entity_count']}",
        "",
        "## Annotation applicability",
        "",
        f"- High-confidence boundary subset: {statistics['high_confidence_intergenic']['count']} / {statistics['intergenic_region_count']} ({statistics['high_confidence_intergenic']['fraction']:.6f})",
        f"- Intervals with at least one uncertain boundary: {statistics['intergenic_region_count'] - statistics['high_confidence_intergenic']['count']} / {statistics['intergenic_region_count']} ({1.0 - statistics['high_confidence_intergenic']['fraction']:.6f})",
        "- Default buffer-distance and minimum-window threshold inference: unavailable",
        "- Scientific acceptance: blocked pending authoritative review of systematic annotation-boundary uncertainty",
        "- These submitted-coordinate intervals remain a descriptive baseline and are not formal candidate windows.",
        "",
        "## Sequence endpoint completeness",
        "",
        "| Sequence | 5′ | 3′ | Evidence |",
        "| --- | --- | --- | --- |",
    ]
    for seqid, values in sorted(statistics["sequence_statistics"].items()):
        endpoints = values["endpoint_completeness"]
        lines.append(
            f"| {seqid} | {endpoints['five_prime']['status']} | {endpoints['three_prime']['status']} | "
            f"5′: {endpoints['five_prime']['evidence_source']}; 3′: {endpoints['three_prime']['evidence_source']} |"
        )
    lines.extend([
        "",
        "## Known assembly and evidence limitations",
        "",
    ])
    lines.extend(f"- {item}" for item in limitations)
    lines.extend(
        [
            "",
            "## Missing risk tracks",
            "",
        ]
    )
    lines.extend(f"- {item}: unavailable" for item in reference["missing_risk_tracks"])
    return "\n".join(lines) + "\n"
