from __future__ import annotations

import json

import streamlit as st

from common import current_entry, load_catalog, render_status_banner
from services.catalog import RUN_TYPE_BASELINE, RUN_TYPE_CANDIDATE_WINDOWS
from services.labels import orientation_label
from services.presentation import load_statistics

st.title("基因组统计 Genome statistics")
st.caption("这里展示的是全基因组基线统计，来自「基因组基线数据」运行本身，不是候选窗口筛选之后的结果。")

entries = load_catalog()
entry = current_entry(entries)
if entry is None:
    st.info("还没有选定运行。请到 **概览 Overview** 页面选择或触发一次运行。")
    st.stop()

baseline_entry = entry
if entry.run_type == RUN_TYPE_CANDIDATE_WINDOWS:
    parent_run_id = entry.manifest.get("parent_run_id")
    baseline_entry = next((item for item in entries if item.run_id == parent_run_id), None)
elif entry.run_type != RUN_TYPE_BASELINE:
    baseline_entry = None

if baseline_entry is None or not baseline_entry.displayable:
    st.error("找不到当前选择对应的可展示基因组基线数据运行。")
    st.stop()

render_status_banner(baseline_entry)

statistics = load_statistics(baseline_entry)

def _relabel_orientation_counts(counts: dict[str, int]) -> dict[str, int]:
    return {orientation_label(key): value for key, value in counts.items()}


columns = st.columns(3)
columns[0].metric("核染色体数量", statistics.get("nuclear_chromosome_count", "n/a"))
columns[1].metric("原始基因间区数量", statistics.get("intergenic_region_count", "n/a"))
high_confidence = statistics.get("high_confidence_intergenic", {})
columns[2].metric(
    "高置信边界子集占比",
    f"{high_confidence.get('count', 'n/a')} ({high_confidence.get('fraction', 0):.2%})"
    if isinstance(high_confidence.get("fraction"), (int, float))
    else "n/a",
)

st.subheader("方向分布（全部提交坐标区间）")
st.bar_chart(_relabel_orientation_counts(statistics.get("orientation_counts", {})))

if high_confidence.get("orientation_counts"):
    st.subheader("方向分布（仅高置信边界子集）")
    st.bar_chart(_relabel_orientation_counts(high_confidence["orientation_counts"]))

st.subheader("基因间区长度分布 Intergenic length distribution")
st.json(statistics.get("intergenic_length_distribution", {}))

with st.expander("完整 statistics.json"):
    st.json(statistics)
