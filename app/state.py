from __future__ import annotations

import streamlit as st

# Single authoritative declaration of every st.session_state key this app uses
# (ARCHITECTURE.md 4.11). Namespaced so future pages cannot collide by accident.
NAMESPACE = "psh"

KEY_CATALOG_ROOT = f"{NAMESPACE}.catalog_root"
KEY_CURRENT_RUN_ID = f"{NAMESPACE}.current_run_id"
KEY_CANDIDATE_FILTERS = f"{NAMESPACE}.candidate_filters"
KEY_SELECTED_CANDIDATE_ID = f"{NAMESPACE}.selected_candidate_id"
KEY_LAST_TRIGGER_RESULT = f"{NAMESPACE}.last_trigger_result"

DEFAULT_CANDIDATE_FILTERS = {"seqid": "(all)", "structural_tier": "(all)", "orientation": "(all)"}


def ensure_defaults() -> None:
    st.session_state.setdefault(KEY_CATALOG_ROOT, None)
    st.session_state.setdefault(KEY_CURRENT_RUN_ID, None)
    st.session_state.setdefault(KEY_CANDIDATE_FILTERS, dict(DEFAULT_CANDIDATE_FILTERS))
    st.session_state.setdefault(KEY_SELECTED_CANDIDATE_ID, None)
    st.session_state.setdefault(KEY_LAST_TRIGGER_RESULT, None)


def set_current_run(run_id: str | None) -> None:
    """Switch the active run and reset all state derived from the previous one.

    ARCHITECTURE.md 4.11 requires clearing the candidate selection, filters,
    and any other view tied to the old run whenever the current run changes --
    otherwise a filter or selection from one run could silently leak into the
    display of a different run.
    """
    if st.session_state.get(KEY_CURRENT_RUN_ID) == run_id:
        return
    st.session_state[KEY_CURRENT_RUN_ID] = run_id
    st.session_state[KEY_CANDIDATE_FILTERS] = dict(DEFAULT_CANDIDATE_FILTERS)
    st.session_state[KEY_SELECTED_CANDIDATE_ID] = None
