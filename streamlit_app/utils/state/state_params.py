"""
Workspace ``params.json`` I/O — load and save widget state per workspace.

Pure I/O only: no Streamlit widgets. The sidebar / page-setup helpers that
surface these values live in :mod:`streamlit_app.utils.state.state_sidebar`.

Tunable defaults are owned by :data:`core.common.defaults.DEFAULTS`.
When a key is missing from ``params.json``, callers should fall back to
``DEFAULTS.<key>`` at the call site (see ``get_param``).
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import streamlit as st

from streamlit_app.utils.state.state_core import _persisted_key


# Canonical list of workspace params persisted in ``params.json``.
# ``save_params`` pulls each of these from session_state (or its
# ``_persisted`` backup) so widget values survive across page
# navigation and browser reloads even on workspaces whose params.json
# does not yet contain the key. Keep this list in sync with the widget
# keys used in the Setup page, sidebar, and file browser.
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
    """Load the workspace's ``params.json`` — or an empty dict if none exists.

    Missing keys are not populated here; callers resolve them via
    ``get_param(params, key, DEFAULTS.<key>)`` so there is exactly one
    source of truth for default values (``core.common.defaults``).
    """
    # Short-circuit: workspaces disabled or caller asked for defaults
    if default or not st.session_state.settings.get("enable_workspaces", False):
        return {}

    # Read: load params.json from the current workspace
    path = Path(st.session_state.workspace, "params.json")
    if not path.exists():
        return {}

    with open(path, "r", encoding="utf-8") as handle:
        raw_params = json.load(handle)

    # Migrate: hyphenated legacy keys -> snake_case (snake_case wins on dupes)
    # Migrate legacy hyphenated top-level keys (e.g. ``reference-substring``)
    # to their snake_case equivalent. Defaults and in-code keys are snake_case
    # now; without this, workspaces created before the migration keep stale
    # hyphen keys and ``save_params`` can't sync widget values back to disk.
    # Snake_case wins if both forms exist.
    migrated: dict[str, Any] = {key: value for key, value in raw_params.items() if "-" not in key}
    for key, value in raw_params.items():
        if "-" in key:
            migrated.setdefault(key.replace("-", "_"), value)
    return migrated


def save_params(params: dict[str, Any]) -> None:
    """Sync widget values from session_state into params and write to workspace.

    For each param key, checks:
    1. The key itself in session_state (widget value, present on the page that owns it)
    2. A persisted non-widget key (``peptide_mirror.<key>_persisted``), which
       survives page navigation even after Streamlit cleans up widget keys.
    """
    # Short-circuit: workspaces disabled means nothing to persist
    if not st.session_state.settings.get("enable_workspaces", False):
        return

    # Sync: pull live widget values (or _persisted backups) into the params dict
    # Union of the canonical param keys and whatever is already in
    # ``params``. A key only appears in ``params`` if load_params picked it
    # up from params.json; canonical keys need to be saved even on their
    # first-ever write, so iterating ``params`` alone would strand fresh
    # widget values.
    tracked_keys = list(dict.fromkeys((*_PARAM_KEYS, *params.keys())))
    for key in tracked_keys:
        if key in st.session_state:
            params[key] = st.session_state[key]
        else:
            # Fall back to the persisted non-widget copy (set by setup.py
            # fragments) so values survive page navigation.
            persisted = _persisted_key(key)
            if persisted in st.session_state:
                params[key] = st.session_state[persisted]

    # Atomic write: stage to a unique tmp file, then os.replace into place
    # Atomic write via unique staging file. A single Streamlit rerun calls
    # save_params from several places (sidebar, each setup fragment, ...);
    # sharing one ``params.json.tmp`` path between them — or between browser
    # tabs on the same workspace — caused the second mover's os.replace to
    # fail because the first had already moved the staging file away.
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
