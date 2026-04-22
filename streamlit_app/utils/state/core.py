"""
Core session-state helpers shared by ``state.workspace`` and ``state.params``.

These are pure utilities that don't touch Streamlit widgets directly — they
just read/write the namespaced keys our app uses to survive navigation.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import streamlit as st

from mspecviz.common.constants import STATE_PREFIX  # noqa: F401

OS_PLATFORM = sys.platform

# Project root = three levels up from streamlit_app/utils/state/core.py.
# All CWD-relative paths (workspaces dir, assets, settings) should be
# anchored here so that the app behaves identically regardless of the
# directory from which ``streamlit run`` is invoked.
PROJECT_ROOT = Path(__file__).resolve().parents[3]


def get_state_key(key: str) -> str:
    """Return a namespaced session state key (e.g. ``peptide_mirror.foo``)."""
    return f"{STATE_PREFIX}.{key}"


def _persisted_key(key: str) -> str:
    """Non-widget backup key that survives page navigation.

    Widget keys (``'ion_types'``) get garbage-collected by Streamlit when the
    owning page unmounts. Fragments mirror the value into a second key with a
    ``_persisted`` suffix; this helper builds that key so the convention lives
    in exactly one place.
    """
    return get_state_key(key + "_persisted")


def get_param(params: dict, key: str, default: Any = None) -> Any:
    """Resolve a param: widget session_state > persisted backup > params.json."""
    if key in st.session_state:
        return st.session_state[key]
    persisted = _persisted_key(key)
    if persisted in st.session_state:
        return st.session_state[persisted]
    return params.get(key, default)
