from __future__ import annotations

import pandas as pd
import streamlit as st

from common import load_catalog, render_status_banner, require_current_run
from services.catalog import RUN_TYPE_CANDIDATE_WINDOWS
from services.labels import rule_id_label
from services.presentation import load_excluded_regions

st.title("排除记录 Excluded regions")

entries = load_catalog()
entry = require_current_run(entries, expected_type=RUN_TYPE_CANDIDATE_WINDOWS)
if entry is None:
    st.stop()

render_status_banner(entry)

try:
    records = load_excluded_regions(entry)
except (OSError, ValueError) as exc:
    st.error(f"无法加载排除记录：{exc}")
    st.stop()

if not records:
    st.info("这次运行没有排除任何区域。")
    st.stop()

rule_ids = ["(全部)"] + sorted({record["rule_id"] for record in records})
chosen_rule = st.selectbox(
    "排除规则",
    rule_ids,
    format_func=lambda value: value if value == "(全部)" else rule_id_label(value),
)
filtered = records if chosen_rule == "(全部)" else [record for record in records if record["rule_id"] == chosen_rule]

st.write(f"显示 {len(filtered)} / {len(records)} 条排除记录")
st.caption("「详细原因」一列由引擎按具体数值动态生成，保留英文原文以便精确核对，未做翻译。")
table = pd.DataFrame(filtered)[["region_id", "parent_region_id", "seqid", "start", "end", "length", "rule_id", "reason"]]
display_table = table.copy()
display_table["rule_id"] = display_table["rule_id"].map(rule_id_label)
st.dataframe(
    display_table,
    width="stretch",
    hide_index=True,
    column_config={
        "region_id": "排除记录编号",
        "parent_region_id": "所属基因间区编号",
        "seqid": "染色体/序列编号",
        "start": st.column_config.NumberColumn("起始位置"),
        "end": st.column_config.NumberColumn("终止位置"),
        "length": st.column_config.NumberColumn("长度（bp）"),
        "rule_id": "排除规则",
        "reason": "详细原因（英文）",
    },
)

st.download_button(
    "下载筛选结果 CSV",
    data=table.to_csv(index=False),
    file_name=f"{entry.run_id}_excluded_regions_filtered.csv",
    mime="text/csv",
)
