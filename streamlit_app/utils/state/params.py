"""
Params I/O and page setup — the stateful entry points used by every page.

``page_setup()`` is called at the top of each view to configure the page,
initialise the workspace, and render the shared sidebar. ``load_params`` /
``save_params`` bridge workspace ``params.json`` with Streamlit session state.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import streamlit as st

from streamlit_app.utils.state.core import PROJECT_ROOT, _persisted_key
from streamlit_app.utils.state.workspace import (
    _get_workspaces_dir,
    _render_workspace_switcher,
    is_safe_workspace_name,
)


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


def page_setup(page: str = "") -> dict[str, Any]:
    """Set up page config, workspace, and sidebar. Call at the top of every page."""
    if "settings" not in st.session_state:
        settings_path = PROJECT_ROOT / "streamlit_app" / "config" / "settings.json"
        with open(settings_path, "r") as f:
            st.session_state.settings = json.load(f)

    logo_path = PROJECT_ROOT / "assets" / "openms_transparent_bg_logo.svg"
    st.set_page_config(
        page_title=st.session_state.settings.get("app-name", "MSpecViz"),
        page_icon=str(logo_path) if logo_path.exists() else None,
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    if "workspace" not in st.session_state or (
        "workspace" in st.query_params
        and st.query_params.workspace != st.session_state.workspace.name
    ):
        st.cache_data.clear()
        st.cache_resource.clear()
        st.session_state.location = "local"

        workspaces_dir = _get_workspaces_dir()
        if st.session_state.settings.get("enable_workspaces", False):
            if "workspace" in st.query_params:
                requested_workspace = st.query_params.workspace
                if not is_safe_workspace_name(requested_workspace):
                    st.session_state.workspace = Path(workspaces_dir, "default")
                    st.query_params.workspace = "default"
                else:
                    st.session_state.workspace = Path(workspaces_dir, requested_workspace)
            else:
                st.session_state.workspace = Path(workspaces_dir, "default")
                st.query_params.workspace = "default"
        else:
            st.session_state.workspace = Path(workspaces_dir, "default")

    if (
        "workspace" not in st.query_params
        and st.session_state.settings.get("enable_workspaces", False)
    ):
        st.query_params.workspace = st.session_state.workspace.name

    st.session_state.workspace.mkdir(parents=True, exist_ok=True)

    return render_sidebar(page)


def render_sidebar(page: str = "") -> dict[str, Any]:
    """Render sidebar with workspace switcher, settings, and version info."""
    params = load_params()

    with st.sidebar:
        if st.session_state.settings.get("enable_workspaces", False):
            _render_workspace_switcher(params)

        with st.expander("**Settings**"):
            img_formats = ["svg", "png", "jpeg", "webp"]
            default_format = params.get("image-format", "png")
            if default_format not in img_formats:
                default_format = "png"
            st.selectbox(
                "image export format",
                img_formats,
                img_formats.index(default_format),
                key="image-format",
            )

        version_info = st.session_state.settings.get("version", "0.1.0")
        app_name = st.session_state.settings.get("app-name", "MSpecViz")
        st.markdown(
            '<style>.version-box{border:1px solid #a4a5ad;padding:10px;'
            'border-radius:.5rem;text-align:center;display:flex;'
            'justify-content:center;align-items:center}</style>',
            unsafe_allow_html=True,
        )
        st.markdown(
            f'<div class="version-box">{app_name}<br>Version: {version_info}</div>',
            unsafe_allow_html=True,
        )

    return params
