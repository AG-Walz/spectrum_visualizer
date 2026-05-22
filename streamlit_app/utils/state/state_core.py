"""Namespaced st.session_state key helpers — imported by every state submodule, view, and page_state."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import streamlit as st

# Prefix on every app-owned session_state key, so they cannot collide with third-party widget keys.
STATE_PREFIX = "peptide_mirror"

# Anchor for CWD-relative paths (workspaces, assets, settings) so the app behaves the
# same no matter which directory streamlit run was launched from.
PROJECT_ROOT = Path(__file__).resolve().parents[3]


def get_state_key(key: str) -> str:
    """Prefix key with the app namespace (e.g. "foo" -> "peptide_mirror.foo")."""
    return f"{STATE_PREFIX}.{key}"


def _persisted_key(key: str) -> str:
    """Return the <key>_persisted backup-key name — survives Streamlit's widget garbage collection on page navigation."""
    return get_state_key(key + "_persisted")


def get_param(params: dict, key: str, default: Any = None) -> Any:
    """Resolve a param with fallback order: live widget > _persisted backup > params.json > default."""
    if key in st.session_state:
        return st.session_state[key]
    persisted = _persisted_key(key)
    if persisted in st.session_state:
        return st.session_state[persisted]
    return params.get(key, default)
