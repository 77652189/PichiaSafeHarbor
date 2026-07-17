from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from pichia_safe_harbor.errors import ContractError
from pichia_safe_harbor.transcript_qualification import qualify_transcript_source


def _identity(path: Path) -> dict:
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _source(tmp_path: Path) -> dict:
    archive = tmp_path / "run.sralite"
    archive.write_bytes(b"archive")
    xml = tmp_path / "run.xml"
    xml.write_text(
        "<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>"
        "<EXPERIMENT accession='SRX1'><STUDY_REF accession='SRP1'/><DESIGN><DESIGN_DESCRIPTION>Truseq mRNA stranded HT; RNA Integrity Number greater than 7</DESIGN_DESCRIPTION><SAMPLE_DESCRIPTOR accession='SRS1'/><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>unspecified</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN><PLATFORM><ILLUMINA><INSTRUMENT_MODEL>NextSeq 550</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>"
        "<STUDY accession='SRP1'><IDENTIFIERS><EXTERNAL_ID namespace='BioProject'>PRJNA1</EXTERNAL_ID></IDENTIFIERS><STUDY_LINKS><STUDY_LINK><XREF_LINK><DB>pubmed</DB><ID>1</ID></XREF_LINK></STUDY_LINK></STUDY_LINKS></STUDY>"
        "<SAMPLE accession='SRS1'><IDENTIFIERS><EXTERNAL_ID namespace='BioSample'>SAMN1</EXTERNAL_ID></IDENTIFIERS><SAMPLE_NAME><SCIENTIFIC_NAME>Komagataella phaffii</SCIENTIFIC_NAME></SAMPLE_NAME><DESCRIPTION>Strain-B low density 6 h methanol</DESCRIPTION><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>strain</TAG><VALUE>Strain-B</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE>"
        "<RUN_SET><RUN accession='SRR1' total_spots='2' total_bases='304'><Statistics><Read index='0' average='76'/><Read index='1' average='76'/></Statistics></RUN></RUN_SET>"
        "</EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>",
        encoding="utf-8",
    )
    policy = tmp_path / "policy.html"
    policy.write_text("policy", encoding="utf-8")
    publication = tmp_path / "publication.xml"
    publication.write_text(
        "<article><license>https://creativecommons.org/licenses/by/4.0/</license><body>Strain-B: ATCC 20864; poly(A)-enriched, strand-specific cDNA from triplicate batch cultivations. de novo assembled transcript models were used for initial annotation.</body></article>",
        encoding="utf-8",
    )
    assembly = tmp_path / "assembly.json"
    assembly.write_text(
        json.dumps(
            {
                "result": {
                    "uids": ["1"],
                    "1": {
                        "assemblyaccession": "GCA_001746955.1",
                        "assemblyname": "ASM174695v1",
                        "gb_bioprojects": [{"bioprojectaccn": "PRJNA304976"}],
                        "biosource": {"infraspecieslist": [{"sub_type": "strain", "sub_value": "Strain-B"}, {"sub_type": "culture-collection", "sub_value": "ATCC:20864"}]},
                        "submitterorganization": "MIT",
                    },
                }
            }
        ),
        encoding="utf-8",
    )
    archive_identity = _identity(archive)
    archive_md5 = hashlib.md5(archive.read_bytes()).hexdigest()
    archive_identity.update({"archive_md5": archive_md5, "actual_md5": archive_md5})
    return {
        "source_id": "fixture",
        "evidence_type": "strand-specific short-read RNA-seq",
        "study_accession": "SRP1",
        "bioproject_accession": "PRJNA1",
        "experiment_accession": "SRX1",
        "biosample_accession": "SAMN1",
        "sra_sample_accession": "SRS1",
        "run_accession": "SRR1",
        "publication": {"pmid": "1", "pmcid": "PMC1"},
        "organism": "Komagataella phaffii",
        "strain": "Strain-B",
        "strain_identity": {"secondary_identifier_usable": False},
        "experimental_condition": {
            "sample_title": "Strain-B low density 6 h methanol",
            "carbon_source": "methanol",
            "cultivation_time": "6 h",
        },
        "library": {
            "strategy": "RNA-Seq",
            "source": "TRANSCRIPTOMIC",
            "selection": "unspecified",
            "layout": "PAIRED",
            "platform": "ILLUMINA",
            "instrument": "NextSeq 550",
            "read_count_per_spot": 2,
            "mean_read_length_bases": 76,
            "total_spots": 2,
            "total_bases": 304,
            "preparation_kit": "TruSeq mRNA Stranded HT",
            "quality_gate": "RNA Integrity Number greater than 7",
            "strand_specificity": "explicitly stranded by library kit metadata",
            "strand_orientation_convention": "unknown; not stated by the archive metadata and not inferred",
        },
        "license": {"file_license_assignment": "not-explicitly-assigned", "database_usage_terms": "NCBI terms", "database_usage_terms_url": "https://example.invalid"},
        "independence": {"relationship": "provenance overlap", "annotation_provenance_overlap": True},
        "coordinate_identity": {"primary_coordinate_compatibility": "pending explicit probe against GCA_001746955.1", "same_study_primary_assembly": "GCA_001746955.1", "same_study_assembly_bioproject": "PRJNA304976", "annotation_provenance_overlap": True},
        "files": {
            "sra_lite": archive_identity,
            "sra_metadata_xml": _identity(xml),
            "policy": _identity(policy),
            "publication_fulltext": _identity(publication),
            "primary_assembly_summary": _identity(assembly),
        },
    }


def test_transcript_source_file_and_metadata_identity_pass(tmp_path: Path) -> None:
    result = qualify_transcript_source(_source(tmp_path), Path("."))
    assert result["qualification_status"] == "file-identity-qualified_independence-license-mapping-pending"
    assert result["gates"]["primary_coordinate_compatibility"] == "pending"
    assert result["gates"]["license"] == "unresolved"
    assert result["gates"]["independence_from_current_annotation"] == "not-passed"
    assert result["eligible_for_full_genome_boundary_track"] is False


@pytest.mark.parametrize("field", ["strain", "license", "strand"])
def test_transcript_source_rejects_unknown_critical_metadata(tmp_path: Path, field: str) -> None:
    source = _source(tmp_path)
    if field == "strain":
        source["strain"] = "unknown"
    elif field == "license":
        source["license"] = {}
    else:
        source["library"]["strand_specificity"] = "unknown"
    with pytest.raises(ContractError):
        qualify_transcript_source(source, Path("."))


def test_transcript_source_rejects_file_identity_mismatch(tmp_path: Path) -> None:
    source = _source(tmp_path)
    source["files"]["sra_lite"]["sha256"] = "0" * 64
    with pytest.raises(ContractError, match="identity mismatch"):
        qualify_transcript_source(source, Path("."))
