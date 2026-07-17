from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from .errors import ContractError
from .io_utils import sha256_file


def load_transcript_source_manifest(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema_version") != 1 or not value.get("sources"):
        raise ContractError("unsupported or empty transcript source manifest")
    if value.get("primary_coordinate_space") != "GCA_001746955.1":
        raise ContractError("transcript source manifest changed the locked primary coordinate space")
    if value.get("exact_target_strain_coordinates") is not False:
        raise ContractError("transcript source manifest must not claim exact Strain-T coordinates")
    return value


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_transcript_source_files(source: dict[str, Any], repo_root: Path) -> dict[str, dict[str, Any]]:
    verified: dict[str, dict[str, Any]] = {}
    for role, identity in sorted(source.get("files", {}).items()):
        path = repo_root / identity["path"]
        if not path.is_file():
            raise ContractError(f"transcript source file is missing: {source['source_id']}:{role}")
        actual = {
            "path": identity["path"],
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        if actual["size_bytes"] != identity["size_bytes"] or actual["sha256"] != identity["sha256"]:
            raise ContractError(f"transcript source file identity mismatch: {source['source_id']}:{role}")
        if identity.get("actual_md5"):
            actual["md5"] = _md5_file(path)
            if actual["md5"] != identity["actual_md5"] or actual["md5"] != identity.get("archive_md5"):
                raise ContractError(f"transcript archive MD5 mismatch: {source['source_id']}:{role}")
        verified[role] = actual
    if not verified:
        raise ContractError(f"transcript source has no verifiable files: {source['source_id']}")
    return verified


def verify_sra_metadata_contract(source: dict[str, Any], metadata_path: Path) -> dict[str, Any]:
    root = ET.parse(metadata_path).getroot()
    experiment = root.find(".//EXPERIMENT")
    study = root.find(".//STUDY")
    run = root.find(".//RUN")
    sample = root.find(".//SAMPLE")
    sample_descriptor = root.find(".//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR")
    design = root.findtext(".//DESIGN_DESCRIPTION", default="")
    sample_description = root.findtext(".//SAMPLE/DESCRIPTION", default="")
    attributes = {
        item.findtext("TAG", default=""): item.findtext("VALUE", default="")
        for item in root.findall(".//SAMPLE_ATTRIBUTE")
    }
    if experiment is None or experiment.get("accession") != source["experiment_accession"]:
        raise ContractError(f"SRA experiment identity mismatch: {source['source_id']}")
    if study is None or study.get("accession") != source["study_accession"]:
        raise ContractError(f"SRA study identity mismatch: {source['source_id']}")
    if run is None or run.get("accession") != source["run_accession"]:
        raise ContractError(f"SRA run identity mismatch: {source['source_id']}")
    if sample is None or sample.get("accession") != source["sra_sample_accession"] or attributes.get("strain") != source["strain"]:
        raise ContractError(f"SRA strain identity mismatch: {source['source_id']}")
    if sample_descriptor is None or sample_descriptor.get("accession") != source["sra_sample_accession"]:
        raise ContractError(f"SRA sample descriptor mismatch: {source['source_id']}")
    external_ids = {
        (item.get("namespace") or "").lower(): (item.text or "")
        for item in root.findall(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
    }
    study_external_ids = {
        (item.get("namespace") or "").lower(): (item.text or "")
        for item in root.findall(".//STUDY/IDENTIFIERS/EXTERNAL_ID")
    }
    if external_ids.get("biosample") != source["biosample_accession"] or study_external_ids.get("bioproject") != source["bioproject_accession"]:
        raise ContractError(f"SRA external accession mismatch: {source['source_id']}")
    organism = root.findtext(".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME", default="")
    if organism != source["organism"]:
        raise ContractError(f"SRA organism identity mismatch: {source['source_id']}")
    if source["experimental_condition"]["sample_title"] not in sample_description:
        raise ContractError(f"SRA experimental condition mismatch: {source['source_id']}")
    quality_gate = source["library"]["quality_gate"]
    design_lower = design.lower()
    quality_supported = quality_gate.lower() in design_lower or ("rin" in design_lower and "> 7" in design)
    if source["library"]["preparation_kit"].lower() not in design_lower or not quality_supported:
        raise ContractError(f"SRA library protocol mismatch: {source['source_id']}")
    library = source["library"]
    strategy = root.findtext(".//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY", default="")
    library_source = root.findtext(".//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE", default="")
    selection = root.findtext(".//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION", default="")
    layout = "PAIRED" if root.find(".//LIBRARY_LAYOUT/PAIRED") is not None else "SINGLE"
    platform = "ILLUMINA" if root.find(".//PLATFORM/ILLUMINA") is not None else "unknown"
    instrument = root.findtext(".//PLATFORM/ILLUMINA/INSTRUMENT_MODEL", default="")
    if (strategy, library_source, selection, layout, platform, instrument) != (
        library["strategy"],
        library["source"],
        library["selection"],
        library["layout"],
        library["platform"],
        library["instrument"],
    ):
        raise ContractError(f"SRA library descriptor mismatch: {source['source_id']}")
    if int(run.get("total_spots", "-1")) != library["total_spots"] or int(run.get("total_bases", "-1")) != library["total_bases"]:
        raise ContractError(f"SRA run totals mismatch: {source['source_id']}")
    read_lengths = [int(item.get("average", "-1")) for item in run.findall("./Statistics/Read")]
    if len(read_lengths) != library["read_count_per_spot"] or set(read_lengths) != {library["mean_read_length_bases"]}:
        raise ContractError(f"SRA read layout mismatch: {source['source_id']}")
    pubmed_ids = [item.text or "" for item in root.findall(".//STUDY_LINK/XREF_LINK[DB='pubmed']/ID")]
    if source["publication"]["pmid"] not in {value.strip() for value in pubmed_ids}:
        raise ContractError(f"SRA publication link mismatch: {source['source_id']}")
    return {
        "study_accession": study.get("accession"),
        "bioproject_accession": study_external_ids["bioproject"],
        "experiment_accession": experiment.get("accession"),
        "run_accession": run.get("accession"),
        "biosample_accession": external_ids["biosample"],
        "sra_sample_accession": sample.get("accession"),
        "organism": organism,
        "strain": attributes["strain"],
        "sample_description": sample_description,
        "library_protocol": design,
        "library_strategy": strategy,
        "library_source": library_source,
        "library_selection": selection,
        "library_layout": layout,
        "platform": platform,
        "instrument": instrument,
        "read_lengths": read_lengths,
        "total_spots": int(run.get("total_spots", "0")),
        "total_bases": int(run.get("total_bases", "0")),
        "publication_pmid": source["publication"]["pmid"],
    }


def verify_publication_and_assembly_context(
    source: dict[str, Any],
    publication_path: Path,
    assembly_summary_path: Path,
) -> dict[str, Any]:
    publication_root = ET.parse(publication_path).getroot()
    publication_text = " ".join("".join(publication_root.itertext()).split())
    publication_lower = publication_text.lower()
    required_phrases = (
        "strain-b: atcc 20864",
        "poly(a)-enriched, strand-specific cdna",
        "triplicate batch cultivations",
        "de novo assembled transcript models",
        "initial annotation",
    )
    if any(phrase not in publication_lower for phrase in required_phrases):
        raise ContractError(f"publication context is incomplete: {source['source_id']}")
    if "creativecommons.org/licenses/by/4.0" not in publication_path.read_text(encoding="utf-8").lower():
        raise ContractError(f"publication license identity mismatch: {source['source_id']}")
    assembly = json.loads(assembly_summary_path.read_text(encoding="utf-8"))
    uids = assembly.get("result", {}).get("uids", [])
    if len(uids) != 1:
        raise ContractError(f"primary assembly query is ambiguous: {source['source_id']}")
    record = assembly["result"][uids[0]]
    bioprojects = {item.get("bioprojectaccn") for item in record.get("gb_bioprojects", [])}
    infraspecies = {
        item.get("sub_type"): item.get("sub_value")
        for item in record.get("biosource", {}).get("infraspecieslist", [])
    }
    coordinate = source["coordinate_identity"]
    if (
        record.get("assemblyaccession") != coordinate["same_study_primary_assembly"]
        or coordinate["same_study_assembly_bioproject"] not in bioprojects
        or infraspecies.get("strain") != "Strain-B"
        or infraspecies.get("culture-collection") != "ATCC:20864"
    ):
        raise ContractError(f"same-study primary assembly identity mismatch: {source['source_id']}")
    if coordinate.get("annotation_provenance_overlap") is not True or source["independence"].get("annotation_provenance_overlap") is not True:
        raise ContractError(f"annotation provenance overlap was not preserved: {source['source_id']}")
    return {
        "publication_pmid": source["publication"]["pmid"],
        "publication_pmcid": source["publication"]["pmcid"],
        "publication_strain": "Strain-B, ATCC 20864",
        "publication_transcript_design": "poly(A)-enriched strand-specific cDNA from triplicate batch cultivations",
        "transcripts_used_for_initial_annotation": True,
        "assembly_accession": record["assemblyaccession"],
        "assembly_name": record["assemblyname"],
        "assembly_bioproject": coordinate["same_study_assembly_bioproject"],
        "assembly_strain": "Strain-B, ATCC 20864",
        "assembly_submitter": record["submitterorganization"],
    }


def qualify_transcript_source(source: dict[str, Any], repo_root: Path) -> dict[str, Any]:
    strain = source.get("strain", "").strip()
    license_record = source.get("license", {})
    condition = source.get("experimental_condition", {})
    library = source.get("library", {})
    if not strain or strain.lower() == "unknown":
        raise ContractError(f"transcript source strain is unknown: {source['source_id']}")
    if not license_record.get("file_license_assignment") or not license_record.get("database_usage_terms") or not license_record.get("database_usage_terms_url"):
        raise ContractError(f"transcript source license metadata is incomplete: {source['source_id']}")
    if not condition.get("sample_title") or not condition.get("carbon_source") or not condition.get("cultivation_time"):
        raise ContractError(f"transcript source experimental condition is incomplete: {source['source_id']}")
    if "stranded" not in library.get("strand_specificity", "").lower():
        raise ContractError(f"transcript source strand specificity is not explicit: {source['source_id']}")
    if library.get("strand_orientation_convention") != "unknown; not stated by the archive metadata and not inferred":
        raise ContractError(f"transcript source strand orientation was inferred: {source['source_id']}")
    verified_files = verify_transcript_source_files(source, repo_root)
    metadata_role = source["files"].get("sra_metadata_xml")
    if metadata_role is None:
        raise ContractError(f"transcript source SRA metadata is missing: {source['source_id']}")
    metadata = verify_sra_metadata_contract(source, repo_root / metadata_role["path"])
    publication_role = source["files"].get("publication_fulltext")
    assembly_role = source["files"].get("primary_assembly_summary")
    if publication_role is None or assembly_role is None:
        raise ContractError(f"transcript source publication or assembly context is missing: {source['source_id']}")
    provenance = verify_publication_and_assembly_context(
        source,
        repo_root / publication_role["path"],
        repo_root / assembly_role["path"],
    )
    coordinate_status = source["coordinate_identity"].get("primary_coordinate_compatibility")
    if coordinate_status != "pending explicit probe against GCA_001746955.1":
        raise ContractError(f"transcript source coordinate status was overclaimed: {source['source_id']}")
    return {
        "source_id": source["source_id"],
        "evidence_type": source["evidence_type"],
        "strain": strain,
        "experimental_condition": condition,
        "library": library,
        "license": source["license"],
        "independence": source["independence"],
        "verified_files": verified_files,
        "verified_metadata": metadata,
        "verified_provenance": provenance,
        "gates": {
            "file_identity": "passed",
            "metadata_identity": "passed",
            "strain_label": "passed",
            "strain_secondary_identifier": "conflict",
            "experimental_condition": "passed",
            "license": "unresolved",
            "strand_specificity": "passed",
            "strand_orientation_convention": "unavailable",
            "independence_from_current_annotation": "not-passed",
            "primary_coordinate_compatibility": "pending",
        },
        "qualification_status": "file-identity-qualified_independence-license-mapping-pending",
        "eligible_for_full_genome_boundary_track": False,
    }
