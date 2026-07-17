from __future__ import annotations

import pandas as pd
import streamlit as st

import state as app_state
from common import load_catalog, render_status_banner, require_current_run
from services.catalog import RUN_TYPE_CANDIDATE_WINDOWS
from services.labels import (
    collinearity_status_label,
    evidence_level_label,
    orientation_label,
    rule_id_label,
    structural_tier_label,
)
from services.presentation import candidate_view, load_candidate_windows

app_state.ensure_defaults()

st.title("候选位点 Candidates")

entries = load_catalog()
entry = require_current_run(entries, expected_type=RUN_TYPE_CANDIDATE_WINDOWS)
if entry is None:
    st.stop()

render_status_banner(entry)
st.caption(
    f"这里的共线性状态恒为「{collinearity_status_label('unavailable')}」：跨菌株共线性比对"
    "（Strain-B 与 Strain-C 之间）在当前版本中被明确排除、暂不实现（ADR-0011），不是被遗漏了。"
)

try:
    records = [candidate_view(record) for record in load_candidate_windows(entry)]
except (OSError, ValueError) as exc:
    st.error(f"无法加载候选窗口：{exc}")
    st.stop()

if not records:
    st.info("这次运行没有产生任何候选窗口。")
    st.stop()

filters = st.session_state[app_state.KEY_CANDIDATE_FILTERS]
seqids = ["(全部)"] + sorted({record["seqid"] for record in records})
tiers = ["(全部)"] + sorted({record["structural_tier"] for record in records})
orientations = ["(全部)"] + sorted({record["orientation"] for record in records})

columns = st.columns(3)
filters["seqid"] = columns[0].selectbox(
    "染色体/序列",
    seqids,
    index=seqids.index(filters.get("seqid", "(全部)")) if filters.get("seqid") in seqids else 0,
)
filters["structural_tier"] = columns[1].selectbox(
    "结构分层",
    tiers,
    index=tiers.index(filters.get("structural_tier", "(全部)")) if filters.get("structural_tier") in tiers else 0,
    format_func=lambda value: value if value == "(全部)" else structural_tier_label(value),
)
filters["orientation"] = columns[2].selectbox(
    "邻近基因方向",
    orientations,
    index=orientations.index(filters.get("orientation", "(全部)")) if filters.get("orientation") in orientations else 0,
    format_func=lambda value: value if value == "(全部)" else orientation_label(value),
)
st.session_state[app_state.KEY_CANDIDATE_FILTERS] = filters

filtered = [
    record
    for record in records
    if (filters["seqid"] == "(全部)" or record["seqid"] == filters["seqid"])
    and (filters["structural_tier"] == "(全部)" or record["structural_tier"] == filters["structural_tier"])
    and (filters["orientation"] == "(全部)" or record["orientation"] == filters["orientation"])
]

def _translate_rule_flags(raw: str) -> str:
    if not raw:
        return raw
    return "、".join(rule_id_label(rule_id.strip()) for rule_id in raw.split(","))


st.write(f"显示 {len(filtered)} / {len(records)} 个候选窗口")
table = pd.DataFrame(filtered)[
    [
        "candidate_id",
        "seqid",
        "start",
        "end",
        "length",
        "orientation",
        "structural_tier",
        "evidence_level",
        "collinearity_status",
        "rule_flag_ids",
    ]
]
# table keeps raw engine values (used for the CSV export below); display_table
# is a translated copy for the on-screen view only.
display_table = table.copy()
display_table["orientation"] = display_table["orientation"].map(orientation_label)
display_table["structural_tier"] = display_table["structural_tier"].map(structural_tier_label)
display_table["evidence_level"] = display_table["evidence_level"].map(evidence_level_label)
display_table["collinearity_status"] = display_table["collinearity_status"].map(collinearity_status_label)
display_table["rule_flag_ids"] = display_table["rule_flag_ids"].map(_translate_rule_flags)
st.dataframe(
    display_table,
    width="stretch",
    hide_index=True,
    column_config={
        "candidate_id": "候选编号",
        "seqid": "染色体/序列编号",
        "start": st.column_config.NumberColumn("起始位置"),
        "end": st.column_config.NumberColumn("终止位置"),
        "length": st.column_config.NumberColumn("长度（bp）"),
        "orientation": "邻近基因方向",
        "structural_tier": "结构分层",
        "evidence_level": "证据等级",
        "collinearity_status": "共线性状态",
        "rule_flag_ids": "命中规则",
    },
)

candidate_ids = [record["candidate_id"] for record in filtered]
selected = st.session_state.get(app_state.KEY_SELECTED_CANDIDATE_ID)
if selected not in candidate_ids:
    selected = candidate_ids[0] if candidate_ids else None
if candidate_ids:
    selected = st.selectbox("查看某个候选的详情", candidate_ids, index=candidate_ids.index(selected))
    st.session_state[app_state.KEY_SELECTED_CANDIDATE_ID] = selected
    detail = next(record for record in filtered if record["candidate_id"] == selected)

    st.subheader("候选窗口序列")
    sequence = detail.get("sequence", "")
    if sequence:
        st.caption(f"长度 {len(sequence)} bp，可直接复制用于 BLAST 比对、引物/同源臂设计等后续分析。")
        st.code(sequence, language=None)
        st.download_button(
            "下载这条序列（FASTA）",
            data=f">{detail['candidate_id']} seqid={detail['seqid']} start={detail['start']} end={detail['end']} length={detail['length']}\n{sequence}\n",
            file_name=f"{detail['candidate_id'].replace(':', '_')}.fasta",
            mime="text/plain",
            key=f"download-sequence-{detail['candidate_id']}",
        )
    else:
        st.info("这条候选没有记录序列。")

    with st.expander("完整原始数据（英文字段名，供排查问题参考）", expanded=False):
        st.json(detail)

download_columns = st.columns(2)
download_columns[0].download_button(
    "下载筛选结果 CSV",
    data=table.to_csv(index=False),
    file_name=f"{entry.run_id}_candidate_windows_filtered.csv",
    mime="text/csv",
)
fasta_path = entry.run_dir / "candidate_windows.fasta"
if fasta_path.is_file():
    download_columns[1].download_button(
        "下载全部候选序列（FASTA）",
        data=fasta_path.read_bytes(),
        file_name=f"{entry.run_id}_candidate_windows.fasta",
        mime="text/plain",
    )
