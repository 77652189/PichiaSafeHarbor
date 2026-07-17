from __future__ import annotations

import streamlit as st

from common import current_entry, load_catalog, render_status_banner
from services.catalog import RUN_TYPE_BASELINE
from services.labels import collinearity_status_label, missing_risk_track_label
from services.presentation import load_report_text

st.title("数据与方法 Data and methods")

entries = load_catalog()
entry = current_entry(entries)
if entry is None:
    st.info("还没有选定运行。请到 **概览 Overview** 页面选择或触发一次运行。")
    st.stop()

render_status_banner(entry)
manifest = entry.manifest

st.subheader("坐标约定")
st.write(
    f"- 坐标系统：`{manifest.get('coordinate_system', 'n/a')}`"
    "（0 起始、左闭右开区间）"
)
st.write(f"- 主参考基因组版本：`{manifest.get('primary_assembly', 'n/a')}`")
if manifest.get("secondary_assembly"):
    st.write(
        f"- 辅助参考基因组版本：`{manifest.get('secondary_assembly')}`"
        f"（共线性核查：{collinearity_status_label('unavailable')} —— "
        "跨菌株共线性比对功能尚未实现，见 ADR-0011）"
    )

if entry.run_type == RUN_TYPE_BASELINE:
    st.subheader("输入文件")
    for role, identity in manifest.get("inputs", {}).items():
        st.write(f"- `{role}`：`{identity.get('path')}` —— sha256 `{identity.get('sha256')}`（{identity.get('size_bytes')} 字节）")
else:
    st.subheader("规则版本与参数")
    st.write(f"- 规则版本 rule_version：`{manifest.get('rule_version', 'n/a')}`")
    st.json(manifest.get("rule_params", {}))
    st.write(f"- 父运行 parent_run_id：`{manifest.get('parent_run_id', 'n/a')}`")

st.subheader("已知限制")
st.caption("以下条目来自参考基因组数据来源（含 NCBI 官方元数据），保留英文原文以确保准确、可追溯。")
for item in manifest.get("known_limitations", []):
    st.write(f"- {item}")

st.subheader("缺失的风险轨道")
for item in manifest.get("missing_risk_tracks", []):
    st.write(f"- {missing_risk_track_label(item)}：{collinearity_status_label('unavailable')}")

st.subheader("软件版本")
software = manifest.get("software", {})
st.write(f"- {software.get('name', 'n/a')} {software.get('version', 'n/a')}")
st.write(f"- 实现代码哈希 implementation_sha256：`{software.get('implementation_sha256', 'n/a')}`")

report_text = load_report_text(entry)
if report_text:
    with st.expander("完整运行报告 Full run report"):
        st.markdown(report_text)

with st.expander("完整 run_manifest.json"):
    st.json(manifest)
if entry.acceptance is not None:
    with st.expander("完整 acceptance_manifest.json"):
        st.json(entry.acceptance)
