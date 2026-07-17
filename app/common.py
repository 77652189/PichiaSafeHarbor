from __future__ import annotations

from pathlib import Path

import streamlit as st

import state as app_state
from services.catalog import RunCatalogEntry, run_type_label, scan_run_catalog
from services.labels import (
    execution_status_label,
    run_purpose_label,
    scientific_acceptance_blocker_label,
    scientific_acceptance_status_label,
    verification_status_label,
)
from services.presentation import status_banner

DEFAULT_CATALOG_ROOT = Path(__file__).resolve().parent.parent / "local_runs" / "streamlit_runs"


def catalog_root() -> Path:
    app_state.ensure_defaults()
    stored = st.session_state.get(app_state.KEY_CATALOG_ROOT)
    return Path(stored) if stored else DEFAULT_CATALOG_ROOT


def load_catalog() -> list[RunCatalogEntry]:
    """Deliberately uncached: scanning a handful of small run_manifest.json files
    is cheap, and ARCHITECTURE.md 4.11 forbids a stale, unversioned cross-run
    cache far more strongly than it requires caching this. Recomputing every
    rerun means a just-triggered run is always visible immediately."""
    return scan_run_catalog(catalog_root())


def clear_catalog_cache() -> None:
    """No-op now that load_catalog() is uncached; kept so callers don't need to
    change if caching is reintroduced later with a correct invalidation story."""


def current_entry(entries: list[RunCatalogEntry]) -> RunCatalogEntry | None:
    app_state.ensure_defaults()
    run_id = st.session_state.get(app_state.KEY_CURRENT_RUN_ID)
    if not run_id:
        return None
    return next((entry for entry in entries if entry.run_id == run_id), None)


def render_status_banner(entry: RunCatalogEntry) -> None:
    banner = status_banner(entry)
    if banner["scientific_acceptance_status"] != "accepted":
        st.warning(
            f"科学验收状态：**{scientific_acceptance_status_label(banner['scientific_acceptance_status'])}** —— "
            "这次运行不是已验证的科学结论。",
            icon="⚠️",
        )
    columns = st.columns(4)
    columns[0].metric("执行状态", execution_status_label(banner["execution_status"]))
    columns[1].metric("验证状态", verification_status_label(banner["verification_status"]))
    columns[2].metric("科学验收状态", scientific_acceptance_status_label(banner["scientific_acceptance_status"]))
    strain_label = banner["target_strain"]
    if not banner["exact_target_strain_coordinates"]:
        strain_label += "（近缘代理，非精确坐标）"
    columns[3].metric("目标菌株", strain_label)
    if banner.get("run_purpose"):
        st.caption(f"运行用途：{run_purpose_label(banner['run_purpose'])}")
    blockers = banner.get("scientific_acceptance_blockers") or []
    if blockers:
        with st.expander("为什么这次运行没有通过科学验收"):
            for blocker in blockers:
                st.write(f"- {scientific_acceptance_blocker_label(blocker)}")


def require_current_run(entries: list[RunCatalogEntry], expected_type: str | None = None) -> RunCatalogEntry | None:
    entry = current_entry(entries)
    if entry is None:
        st.info("还没有选定运行。请到 **概览 Overview** 页面选择或触发一次运行。")
        return None
    if not entry.displayable:
        st.error(f"运行 `{entry.run_id}` 无法展示：{entry.reason}")
        return None
    if expected_type is not None and entry.run_type != expected_type:
        st.info(
            f"当前选定的运行是「{run_type_label(entry.run_type)}」；这个页面需要"
            f"「{run_type_label(expected_type)}」类型的运行。请到 **概览 Overview** 页面重新选择。"
        )
        return None
    return entry
