from __future__ import annotations

from pathlib import Path
import hashlib

import pytest


@pytest.fixture
def fixture_reference(fixture_paths) -> dict:
    fasta_path, gff_path = fixture_paths

    def file_spec(path: Path) -> dict:
        return {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size_bytes": path.stat().st_size}

    return {
        "key": "fixture",
        "manifest_version": "test-v1",
        "assembly_accession": "FIXTURE_001.1",
        "annotation_accession": "fixture annotation v1",
        "require_gff_assembly_directive": True,
        "sequence_name_map": {},
        "sequence_classes": {
            "chr1": "nuclear_chromosome",
            "chr2": "nuclear_chromosome",
            "chr_empty": "nuclear_chromosome",
            "mito": "mitochondrial",
            "scaffold": "unlocalized_scaffold",
        },
        "sequence_endpoints": {
            "chr1": {
                "five_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
                "three_prime": {"status": "partial", "evidence_source": "fixture", "original_declaration": "partial sequence"},
            },
            "chr2": {
                "five_prime": {"status": "unknown", "evidence_source": "fixture", "original_declaration": "unavailable"},
                "three_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
            },
            "chr_empty": {
                "five_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
                "three_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
            },
            "mito": {
                "five_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
                "three_prime": {"status": "complete", "evidence_source": "fixture", "original_declaration": "complete"},
            },
            "scaffold": {
                "five_prime": {"status": "unknown", "evidence_source": "fixture", "original_declaration": "unavailable"},
                "three_prime": {"status": "unknown", "evidence_source": "fixture", "original_declaration": "unavailable"},
            },
        },
        "known_limitations": ["fixture limitation"],
        "missing_risk_tracks": ["fixture risk track"],
        "files": {
            "fasta": file_spec(fasta_path),
            "annotation": file_spec(gff_path),
        },
    }


@pytest.fixture
def fixture_paths() -> tuple[Path, Path]:
    root = Path(__file__).parents[1]
    return root / "testdata/slice0_fixture.fna", root / "testdata/slice0_fixture.gff3"
