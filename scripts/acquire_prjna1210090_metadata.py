from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


SELECTED = {
    "SRR31989016": {"experiment": "SRX27343955", "sample": "SAMN46238833", "replicate": "3"},
    "SRR31989027": {"experiment": "SRX27343944", "sample": "SAMN46238832", "replicate": "2"},
    "SRR31989028": {"experiment": "SRX27343943", "sample": "SAMN46238831", "replicate": "1"},
}
PROTECTED = {
    "slice0": "local_runs/strain-b_slice0_completion_v6_run1/acceptance_manifest.json",
    "slice0a": "local_runs/slice0a/qualification_v5_run1/acceptance_manifest.json",
    "slice0b": "local_runs/slice0b/qualification_v5_run1/acceptance_manifest.json",
    "prjna604658": "local_runs/independent_transcript_sources/prjna604658/qualification_v6_run1/acceptance_manifest.json",
}


def identity(path: Path, relative: str) -> dict:
    return {"path": relative, "size_bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def fetch(url: str) -> tuple[bytes, dict]:
    request = urllib.request.Request(url, headers={"User-Agent": "PichiaSafeHarbor-metadata-qualification/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        data = response.read()
        return data, {"url": url, "http_status": response.status, "content_type": response.headers.get("Content-Type")}


def write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    candidate_dir = root / "transcript_source_candidates/prjna1210090"
    source_dir = root / "local_runs/independent_transcript_sources/prjna1210090/source_files"
    if candidate_dir.exists() or source_dir.exists():
        raise SystemExit("PRJNA1210090 acquisition outputs already exist")
    temp_root = Path(tempfile.mkdtemp(prefix="prjna1210090.", dir=root / "local_runs"))
    try:
        esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "sra", "term": "SRP557139[All Fields]", "retmax": 100, "retmode": "json"})
        esearch_data, esearch_meta = fetch(esearch_url)
        ids = json.loads(esearch_data)["esearchresult"]["idlist"]
        sra_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db": "sra", "id": ",".join(ids), "retmode": "xml"})
        ena_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
            "accession": "SRP557139", "result": "read_run", "format": "json",
            "fields": "study_accession,secondary_study_accession,experiment_accession,sample_accession,secondary_sample_accession,run_accession,library_strategy,library_source,library_selection,library_layout,instrument_model,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp",
        })
        bioproject_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db": "bioproject", "id": "1210090", "retmode": "xml"})
        policy_urls = {
            "ncbi_policy": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
            "ena_terms": "https://www.ebi.ac.uk/about/terms-of-use",
        }
        requests = {"sra_experiment_packages": sra_url, "ena_read_run": ena_url, "bioproject": bioproject_url, **policy_urls}
        snapshots = {}
        request_meta = {"sra_esearch": esearch_meta}
        names = {
            "sra_experiment_packages": "sra_experiment_packages.xml",
            "ena_read_run": "ena_read_run.json",
            "bioproject": "bioproject_1210090.xml",
            "ncbi_policy": "ncbi_policies.html",
            "ena_terms": "ena_terms_of_use.html",
        }
        for key, url in requests.items():
            data, meta = fetch(url)
            path = temp_root / names[key]
            write(path, data)
            snapshots[key] = path
            request_meta[key] = meta
        acquired_at = datetime.now(timezone.utc).isoformat()
        final_source_rel = "local_runs/independent_transcript_sources/prjna1210090/source_files"
        file_identities = {
            key: identity(path, f"{final_source_rel}/{path.name}") for key, path in snapshots.items()
        }
        candidate_dir.mkdir(parents=True)
        acquisition_path = candidate_dir / "acquisition_evidence.v1.json"
        acquisition = {
            "schema_version": 1,
            "candidate_id": "prjna1210090-srp557139-wt-no-stress",
            "acquired_at_utc": acquired_at,
            "requests": request_meta,
            "selected_runs": SELECTED,
            "raw_read_acquisition": {
                "status": "not-performed",
                "reason": "metadata hard gates must pass before downloading approximately 14.8 GB of selected paired FASTQ files",
                "declared_identities_source": "ena_read_run",
            },
            "snapshots": file_identities,
        }
        acquisition_path.write_text(json.dumps(acquisition, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        protected = {}
        for key, relative in PROTECTED.items():
            path = root / relative
            protected[key] = identity(path, relative)
        manifest = {
            "schema_version": 1,
            "candidate_id": "prjna1210090-srp557139-wt-no-stress",
            "accessions": {
                "bioproject": "PRJNA1210090", "study": "SRP557139",
                "runs": SELECTED,
            },
            "reference_annotation_release": "2016-09-21",
            "protected_authority": protected,
            "files": {
                "acquisition_evidence": identity(acquisition_path, "transcript_source_candidates/prjna1210090/acquisition_evidence.v1.json"),
                "reference_manifest": identity(root / "reference/manifest.v1.json", "reference/manifest.v1.json"),
                **file_identities,
            },
        }
        (candidate_dir / "manifest.v1.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        source_dir.parent.mkdir(parents=True, exist_ok=True)
        temp_root.replace(source_dir)
        print(json.dumps({"candidate_id": manifest["candidate_id"], "selected_runs": sorted(SELECTED), "raw_reads": "not-downloaded"}, indent=2))
        return 0
    except Exception:
        if temp_root.exists():
            import shutil
            shutil.rmtree(temp_root, ignore_errors=True)
        if candidate_dir.exists():
            import shutil
            shutil.rmtree(candidate_dir, ignore_errors=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
