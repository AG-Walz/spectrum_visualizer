"""
Workspace ``params.json`` I/O — load and save widget state per workspace.

Pure I/O only: no Streamlit widgets. The sidebar / page-setup helpers that
surface these values live in :mod:`streamlit_app.utils.state.sidebar`.

Tunable defaults are owned by :data:`mspecviz.common.defaults.DEFAULTS`.
When a key is missing from ``params.json``, callers should fall back to
``DEFAULTS.<key>`` at the call site (see ``get_param``).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import streamlit as st

from streamlit_app.utils.state.core import _persisted_key


def load_params(default: bool = False) -> dict[str, Any]:
    """Load the workspace's ``params.json`` — or an empty dict if none exists.

    Missing keys are not populated here; callers resolve them via
    ``get_param(params, key, DEFAULTS.<key>)`` so there is exactly one
    source of truth for default values (``mspecviz.common.defaults``).
    """
    if default or not st.session_state.settings.get("enable_workspaces", False):
        return {}

    path = Path(st.session_state.workspace, "params.json")
    if not path.exists():
        return {}

    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    # Migrate legacy hyphenated top-level keys (e.g. ``reference-substring``)
    # to their snake_case equivalent. Defaults and in-code keys are snake_case
    # now; without this, workspaces created before the migration keep stale
    # hyphen keys and ``save_params`` can't sync widget values back to disk.
    # Snake_case wins if both forms exist.
    migrated: dict[str, Any] = {k: v for k, v in raw.items() if "-" not in k}
    for k, v in raw.items():
        if "-" in k:
            migrated.setdefault(k.replace("-", "_"), v)
    return migrated


def save_params(params: dict[str, Any]) -> None:
    """Sync widget values from session_state into params and write to workspace.

    For each param key, checks:
    1. The key itself in session_state (widget value, present on the page that owns it)
    2. A persisted non-widget key (``peptide_mirror.<key>_persisted``), which
       survives page navigation even after Streamlit cleans up widget keys.
    """
    if not st.session_state.settings.get("enable_workspaces", False):
        return

    for key in params:
        if key in st.session_state:
            params[key] = st.session_state[key]
        else:
            # Fall back to the persisted non-widget copy (set by setup.py
            # fragments) so values survive page navigation.
            persisted = _persisted_key(key)
            if persisted in st.session_state:
                params[key] = st.session_state[persisted]

    # Atomic write: a Streamlit rerun mid-write must never leave a 0-byte file
    # that breaks load_params() on the next page load.
    path = Path(st.session_state.workspace, "params.json")
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w", encoding="utf-8") as outfile:
        json.dump(params, outfile, indent=4)
    os.replace(tmp, path)
