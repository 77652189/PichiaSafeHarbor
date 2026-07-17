from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


NUCLEAR = {"CP014715.1", "CP014716.1", "CP014717.1", "CP014718.1"}
BIN_SIZE = 1000


def median_from_histogram(histogram: Counter, total: int) -> float:
    if total == 0:
        return 0.0
    targets = {(total - 1) // 2, total // 2}
    seen = 0; values = []
    for depth, count in sorted(histogram.items()):
        for target in sorted(targets):
            if seen <= target < seen + count:
                values.append(depth)
        seen += count
    return sum(values) / len(values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows-json", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    parser.add_argument("--windows-tsv", type=Path, required=True)
    parser.add_argument("--anomalies-tsv", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.windows_json.read_text(encoding="utf-8"))
    lengths = config["sequence_lengths"]
    windows = config["windows"]
    stats = {name: {"sum": 0, "covered": 0, "max": 0, "hist": Counter(), "bins": defaultdict(int)} for name in lengths}
    window_acc = {item["window_id"]: {"sum": 0, "covered": 0, "max": 0} for item in windows}
    windows_by_ref = defaultdict(list)
    for item in windows: windows_by_ref[item["reference"]].append(item)
    for line in sys.stdin:
        reference, pos_text, depth_text = line.rstrip("\n").split("\t")[:3]
        if reference not in stats: continue
        pos0 = int(pos_text) - 1; depth = int(depth_text); state = stats[reference]
        state["sum"] += depth; state["covered"] += depth > 0; state["max"] = max(state["max"], depth); state["hist"][depth] += 1
        state["bins"][pos0 // BIN_SIZE] += depth
        for window in windows_by_ref[reference]:
            if window["start_0based"] <= pos0 < window["end_0based"]:
                value = window_acc[window["window_id"]]; value["sum"] += depth; value["covered"] += depth > 0; value["max"] = max(value["max"], depth)
    per_sequence = {}
    anomaly_rows = []
    for reference, length in lengths.items():
        state = stats[reference]
        if sum(state["hist"].values()) != length:
            raise SystemExit(f"depth did not cover every reference position: {reference}")
        bin_means = []
        for index in range((length + BIN_SIZE - 1) // BIN_SIZE):
            width = min(BIN_SIZE, length - index * BIN_SIZE)
            bin_means.append(state["bins"].get(index, 0) / width)
        nonzero = [value for value in bin_means if value > 0]
        median_nonzero = statistics.median(nonzero) if nonzero else 0.0
        mad = statistics.median(abs(value - median_nonzero) for value in nonzero) if nonzero else 0.0
        threshold = median_nonzero + max(10 * mad, 5 * median_nonzero, 10.0)
        for index, value in enumerate(bin_means):
            if value > threshold:
                anomaly_rows.append((reference, index * BIN_SIZE, min(length, (index + 1) * BIN_SIZE), value, median_nonzero, mad, "high_depth_outlier"))
        per_sequence[reference] = {
            "class": "nuclear_chromosome" if reference in NUCLEAR else "non_nuclear_reference",
            "length": length, "covered_bases": state["covered"], "coverage_fraction": state["covered"] / length,
            "mean_depth": state["sum"] / length, "median_depth": median_from_histogram(state["hist"], length), "max_depth": state["max"],
            "nonzero_1kb_bins": len(nonzero), "total_1kb_bins": len(bin_means), "median_nonzero_bin_depth": median_nonzero, "mad_nonzero_bin_depth": mad,
        }
    window_rows = []
    for window in windows:
        value = window_acc[window["window_id"]]; width = window["end_0based"] - window["start_0based"]
        window_rows.append({**window, "covered_bases": value["covered"], "coverage_fraction": value["covered"] / width, "mean_depth": value["sum"] / width, "max_depth": value["max"]})
    summary = {
        "schema_version": 1, "coordinate_space": "Strain-B GCA_001746955.1", "depth_semantics": "samtools depth -aa -d 0 on coordinate-sorted BAM",
        "strand_support": "unavailable", "bin_size": BIN_SIZE, "per_sequence": per_sequence,
        "nuclear_chromosomes_nonzero": all(per_sequence[name]["covered_bases"] > 0 for name in NUCLEAR),
        "fixed_windows_nonzero": sum(item["covered_bases"] > 0 for item in window_rows), "fixed_windows_total": len(window_rows),
        "high_depth_outlier_bins": len(anomaly_rows),
        "nuclear_high_depth_outlier_bins": sum(row[0] in NUCLEAR for row in anomaly_rows),
        "engineering_interpretation": "high RNA-seq depth is expression evidence, not proof of integration; flagged bins are excluded from strain-specific inference pending orthogonal sequence evidence",
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with args.windows_tsv.open("w", encoding="utf-8", newline="\n") as handle:
        fields = ["window_id", "reference", "start_0based", "end_0based", "position_class", "endpoint_caveat", "covered_bases", "coverage_fraction", "mean_depth", "max_depth"]
        handle.write("\t".join(fields) + "\n")
        for item in window_rows: handle.write("\t".join(str(item[field]) for field in fields) + "\n")
    with args.anomalies_tsv.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("reference\tstart_0based\tend_0based\tmean_depth\tchromosome_median_nonzero_bin_depth\tchromosome_mad_nonzero_bin_depth\tclassification\n")
        for row in sorted(anomaly_rows, key=lambda item: (-item[3], item[0], item[1])): handle.write("\t".join(str(value) for value in row) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
