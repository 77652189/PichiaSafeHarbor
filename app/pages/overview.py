from __future__ import annotations

import streamlit as st

import state as app_state
from common import DEFAULT_CATALOG_ROOT, catalog_root, load_catalog, render_status_banner
from services.catalog import run_type_label
from services.labels import execution_status_label, scientific_acceptance_status_label

app_state.ensure_defaults()

st.title("PichiaSafeHarbor —— 概览 Overview")
st.caption(
    "目标菌株：Strain-T。主坐标参考：Strain-B（开发/测试占位数据，非 Strain-T 精确坐标，见 ADR-0009）。"
    "本项目目前没有任何运行通过科学验收。"
)

with st.expander("运行目录扫描位置 Run catalog location", expanded=False):
    current_root = str(catalog_root())
    new_root = st.text_input("扫描运行结果的目录 Directory scanned for runs", value=current_root)
    if new_root != current_root:
        st.session_state[app_state.KEY_CATALOG_ROOT] = new_root
        st.rerun()
    st.caption(f"默认值 Default：`{DEFAULT_CATALOG_ROOT}`")

entries = load_catalog()

if not entries:
    st.info(
        "还没有任何运行结果。请到 **数据导入与运行触发** 页面，先生成一次基因组基线数据，"
        "再用它生成候选安全港窗口。"
    )
else:
    st.caption(
        "不同目录的运行有可能显示相同的 run_id——如果输入完全相同，run_id 是按内容计算的，"
        "内容一样就会一样，这不是错误；下面按运行结果所在的目录名区分。"
    )
    st.subheader("已有运行 Available runs")
    header_cols = st.columns([2, 2, 2, 2, 2, 1])
    for col, label in zip(header_cols, ["运行目录", "类型", "执行状态", "科学验收状态", "run_id", ""]):
        col.markdown(f"**{label}**")
    for entry in entries:
        cols = st.columns([2, 2, 2, 2, 2, 1])
        cols[0].write(f"`{entry.run_dir.name}`")
        cols[1].write(run_type_label(entry.run_type))
        cols[2].write(execution_status_label(entry.execution_status))
        cols[3].write(scientific_acceptance_status_label(entry.scientific_acceptance_status))
        cols[4].write(f"`{entry.run_id}`")
        if entry.displayable:
            if cols[5].button("查看 View", key=f"select-{entry.run_dir.name}"):
                app_state.set_current_run(entry.run_id)
                st.rerun()
        else:
            cols[5].caption(entry.reason)

    current = next(
        (entry for entry in entries if entry.run_id == st.session_state.get(app_state.KEY_CURRENT_RUN_ID)),
        None,
    )
    if current is not None:
        st.subheader(f"当前运行 Current run：`{current.run_dir.name}`（run_id = `{current.run_id}`）")
        render_status_banner(current)
