from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .errors import ContractError
from .acceptance import create_acceptance_manifest
from .pipeline import run_baseline
from .reference import (
    fetch_reference,
    install_reference_archive,
    load_reference_manifest,
    reference_by_key,
    verify_reference_files,
)
from .slice0a import create_slice0a_acceptance, run_slice0a
from .slice1 import create_slice1_acceptance, run_slice1


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pichia-safe-harbor")
    parser.add_argument("--manifest", type=Path, default=Path("reference/manifest.v1.json"))
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name in ("fetch", "validate", "install"):
        command = subparsers.add_parser(name)
        command.add_argument("reference", choices=("strain-b", "strain-c"))
        command.add_argument("--data-dir", type=Path, default=Path("reference/data"))
        if name == "install":
            command.add_argument("--archive", type=Path, required=True)

    baseline = subparsers.add_parser("baseline")
    baseline.add_argument("--reference", choices=("strain-b",), default="strain-b")
    baseline.add_argument("--data-dir", type=Path, default=Path("reference/data"))
    baseline.add_argument("--output-dir", type=Path, required=True)
    acceptance = subparsers.add_parser("acceptance")
    acceptance.add_argument("--run-dir", type=Path, required=True)
    acceptance.add_argument("--independent-check", type=Path, required=True)
    acceptance.add_argument("--test-evidence", type=Path, required=True)
    qualify = subparsers.add_parser("qualify-annotations")
    qualify.add_argument("--sources", type=Path, default=Path("annotation_sources/manifest.v1.json"))
    qualify.add_argument("--output-dir", type=Path, required=True)
    qualify.add_argument("--repo-root", type=Path, default=Path("."))
    qualify_acceptance = subparsers.add_parser("accept-annotation-qualification")
    qualify_acceptance.add_argument("--run-dir", type=Path, required=True)
    qualify_acceptance.add_argument("--independent-check", type=Path, required=True)
    qualify_acceptance.add_argument("--test-evidence", type=Path, required=True)
    candidate_windows = subparsers.add_parser("candidate-windows")
    candidate_windows.add_argument("--baseline-run-dir", type=Path, required=True)
    candidate_windows.add_argument("--reference", choices=("strain-b",), default="strain-b")
    candidate_windows.add_argument("--data-dir", type=Path, default=Path("reference/data"))
    candidate_windows.add_argument("--buffer-distance-bp", type=int, required=True)
    candidate_windows.add_argument("--min-window-bp", type=int, required=True)
    candidate_windows.add_argument("--long-interval-percentile", type=float, default=0.95)
    candidate_windows.add_argument("--output-dir", type=Path, required=True)
    accept_candidate_windows = subparsers.add_parser("accept-candidate-windows")
    accept_candidate_windows.add_argument("--run-dir", type=Path, required=True)
    accept_candidate_windows.add_argument("--independent-check", type=Path, required=True)
    accept_candidate_windows.add_argument("--test-evidence", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "acceptance":
            result = create_acceptance_manifest(
                args.run_dir, args.independent_check, args.test_evidence
            )
            print(
                json.dumps(
                    {
                        "run_id": result["run_id"],
                        "execution_status": result["execution_status"],
                        "verification_status": result["verification_status"],
                        "scientific_acceptance_status": result[
                            "scientific_acceptance_status"
                        ],
                    },
                    indent=2,
                )
            )
            return 0
        if args.command == "qualify-annotations":
            result = run_slice0a(args.sources, args.output_dir, args.repo_root.resolve())
            print(json.dumps({"execution_status": result["execution_status"], "run_id": result["run_id"], "output_dir": str(args.output_dir)}, indent=2))
            return 0
        if args.command == "accept-annotation-qualification":
            result = create_slice0a_acceptance(args.run_dir, args.independent_check, args.test_evidence)
            print(json.dumps({"run_id": result["run_id"], "execution_status": result["execution_status"], "verification_status": result["verification_status"], "scientific_acceptance_status": result["scientific_acceptance_status"]}, indent=2))
            return 0
        if args.command == "accept-candidate-windows":
            result = create_slice1_acceptance(
                args.run_dir, args.independent_check, args.test_evidence
            )
            print(
                json.dumps(
                    {
                        "run_id": result["run_id"],
                        "execution_status": result["execution_status"],
                        "verification_status": result["verification_status"],
                        "scientific_acceptance_status": result["scientific_acceptance_status"],
                    },
                    indent=2,
                )
            )
            return 0
        if args.command == "candidate-windows":
            manifest = load_reference_manifest(args.manifest)
            reference = reference_by_key(manifest, args.reference)
            paths = verify_reference_files(reference, args.data_dir)
            result = run_slice1(
                args.baseline_run_dir,
                reference,
                paths["fasta"],
                args.buffer_distance_bp,
                args.min_window_bp,
                args.output_dir,
                args.long_interval_percentile,
            )
            print(
                json.dumps(
                    {
                        "run_id": result["run_id"],
                        "execution_status": result["execution_status"],
                        "verification_status": result["verification_status"],
                        "scientific_acceptance_status": result["scientific_acceptance_status"],
                        "candidate_window_count": result["candidate_window_count"],
                        "excluded_region_count": result["excluded_region_count"],
                    },
                    indent=2,
                )
            )
            return 0
        manifest = load_reference_manifest(args.manifest)
        reference = reference_by_key(manifest, args.reference)
        if args.command == "fetch":
            paths = fetch_reference(reference, args.data_dir)
            print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2))
        elif args.command == "install":
            paths = install_reference_archive(reference, args.data_dir, args.archive)
            print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2))
        elif args.command == "validate":
            paths = verify_reference_files(reference, args.data_dir)
            print(json.dumps({"status": "valid", "files": {key: str(value) for key, value in paths.items()}}, indent=2))
        else:
            paths = verify_reference_files(reference, args.data_dir)
            result = run_baseline(reference, paths["fasta"], paths["annotation"], args.output_dir)
            print(json.dumps({"execution_status": result["execution_status"], "run_id": result["run_id"], "output_dir": str(args.output_dir)}, indent=2))
        return 0
    except (ContractError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
