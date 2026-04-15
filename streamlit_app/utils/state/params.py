"""
Workspace ``params.json`` I/O — load defaults, save widget state.

Pure I/O only: no Streamlit widgets. The sidebar / page-setup helpers that
surface these values live in :mod:`streamlit_app.utils.state.sidebar`.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import streamlit as st

from streamlit_app.utils.state.core import PROJECT_ROOT, _persisted_key


def load_params(default: bool = False) -> dict[str, Any]:
    """Load default-parameters.json, then overlay workspace params.json on top."""
    default_path = PROJECT_ROOT / "streamlit_app" / "config" / "default-parameters.json"
    with open(default_path, "r", encoding="utf-8") as f:
        params = json.load(f)

    if default or not st.session_state.settings.get("enable_workspaces", False):
        return params

    path = Path(st.session_state.workspace, "params.json")
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            saved = json.load(f)

        # Migrate old format: peptide_mirror.X_y -> x-y
        migrated = False
        if "peptide_mirror" in saved and isinstance(saved["peptide_mirror"], dict):
            old = saved["peptide_mirror"]
            for old_key, value in old.items():
                short = old_key.replace("peptide_mirror.", "")
                new_key = short.replace("_", "-")
                # Prefer old-format value when the flat key is still at its default
                if new_key in params and saved.get(new_key) == params[new_key]:
                    saved[new_key] = value
            # Remove the old dict once migrated
            del saved["peptide_mirror"]
            migrated = True

        # Persist migration so it only runs once
        if migrated:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(saved, f, indent=4)

        # Overlay saved values onto defaults
        for key in params:
            if key in saved:
                params[key] = saved[key]
        # Also restore any extra keys the workspace had (e.g. browser dirs)
        for key in saved:
            if key not in params:
                params[key] = saved[key]

    return params


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
