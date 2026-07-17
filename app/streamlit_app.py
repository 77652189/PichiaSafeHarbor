from __future__ import annotations

import sys
from pathlib import Path

_APP_ROOT = Path(__file__).resolve().parent
_SRC = _APP_ROOT.parent / "src"
for path in (_SRC, _APP_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import streamlit as st  # noqa: E402

st.set_page_config(page_title="PichiaSafeHarbor 安全港候选浏览", layout="wide")

# Single authoritative page/navigation registry (ARCHITECTURE.md 4.10). No other
# module may register a page or duplicate this list.
pages = {
    "概览 Overview": [st.Page("pages/overview.py", title="概览 Overview", icon="🧬")],
    "数据 Data": [
        st.Page("pages/data_import_and_trigger.py", title="数据导入与运行触发", icon="📥"),
        st.Page("pages/data_and_methods.py", title="数据与方法", icon="📋"),
    ],
    "分析 Analysis": [
        st.Page("pages/candidates.py", title="候选位点 Candidates", icon="🧭"),
        st.Page("pages/excluded_regions.py", title="排除记录", icon="🚫"),
        st.Page("pages/genome_statistics.py", title="基因组统计", icon="📊"),
    ],
}
navigation = st.navigation(pages)
navigation.run()
