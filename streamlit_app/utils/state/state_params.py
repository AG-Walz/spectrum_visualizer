"""Load and save each workspace's params.json — called by state_sidebar, page_state, and Setup."""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import streamlit as st

from streamlit_app.utils.state.state_core import _persisted_key


# Canonical list of widget keys save_params writes to params.json. Listing every key
# explicitly (not just those already on disk) lets first-ever saves capture widget values, and
# lets _flush_and_switch_workspace pop them when switching workspaces.
_PARAM_KEYS: tuple[str, ...] = (
    "selected_idxmls", "selected_mzmls",
    "idxml_browser_dir", "mzml_browser_dir",
    "target_peptides", "reference_substring",
    "ion_types", "mz_tolerance",
    "add_neutral_losses",
    "intensity_mode", "workers",
    "koina_enabled", "koina_model",
    "koina_collision_energy", "koina_instrument_type",
    "image_format",
)


def load_params(default: bool = False) -> dict[str, Any]:
    """Load the current workspace's params.json — empty dict if workspaces disabled or no file yet."""
    # Short-circuit: workspaces disabled or caller asked for defaults
    if default or not st.session_state.settings.get("enable_workspaces", False):
        return {}

    # Read: load params.json from the current workspace
    path = Path(st.session_state.workspace, "params.json")
    if not path.exists():
        return {}

    with open(path, "r", encoding="utf-8") as handle:
        raw_params = json.load(handle)

    # Migrate legacy hyphenated keys (reference-substring) to snake_case. Pre-migration
    # workspaces would otherwise hold stale hyphen keys that save_params cannot sync.
    migrated: dict[str, Any] = {key: value for key, value in raw_params.items() if "-" not in key}
    for key, value in raw_params.items():
        if "-" in key:
            migrated.setdefault(key.replace("-", "_"), value)
    return migrated


def save_params(params: dict[str, Any]) -> None:
    """Pull live widget values (or _persisted backups) into params and write atomically to params.json."""
    # Short-circuit: workspaces disabled means nothing to persist
    if not st.session_state.settings.get("enable_workspaces", False):
        return

    # Sync: pull live widget values (or _persisted backups) into the params dict.
    # Union with _PARAM_KEYS so first-ever saves capture widget values not yet on disk.
    tracked_keys = list(dict.fromkeys((*_PARAM_KEYS, *params.keys())))
    for key in tracked_keys:
        if key in st.session_state:
            params[key] = st.session_state[key]
        else:
            # _persisted backup (set by setup.py fragments) — survives widget garbage collection on page navigation.
            persisted = _persisted_key(key)
            if persisted in st.session_state:
                params[key] = st.session_state[persisted]

    # Atomic write via unique staging file. A single rerun can call save_params from the sidebar
    # plus several setup fragments; sharing one params.json.tmp path between them — or between
    # browser tabs on the same workspace — made the second mover's os.replace fail.
    params_path = Path(st.session_state.workspace, "params.json")
    params_path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, staging_name = tempfile.mkstemp(
        prefix=".params-", suffix=".json.tmp", dir=params_path.parent,
    )
    staging_path = Path(staging_name)
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as outfile:
            json.dump(params, outfile, indent=4)
        os.replace(staging_path, params_path)
    except Exception:
        if staging_path.exists():
            try:
                staging_path.unlink()
            except OSError:
                pass
        raise
