"""
Streamlit state and page-setup helpers for MSpecViz.

Adapted from the OpenMS streamlit-template. Provides page_setup(),
render_sidebar(), workspace management, and UI helpers.
"""

import os
import sys
import json
import shutil

import streamlit as st

from typing import Any
from pathlib import Path


OS_PLATFORM = sys.platform


def is_safe_workspace_name(name: str) -> bool:
    """Check if a workspace name is safe (no path traversal characters)."""
    if not name:
        return False
    return "/" not in name and "\\" not in name and name not in ("..", ".")


def load_params(default: bool = False) -> dict[str, Any]:
    """Load default-parameters.json, then overlay workspace params.json on top."""
    default_path = Path(__file__).parents[2] / "default-parameters.json"
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
            # Check for a persisted non-widget copy (set by setup.py fragments)
            persisted_key = f"peptide_mirror.{key.replace('-', '_')}_persisted"
            if persisted_key in st.session_state:
                params[key] = st.session_state[persisted_key]

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
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

    st.set_page_config(
        page_title=st.session_state.settings.get("app-name", "MSpecViz"),
        page_icon="assets/openms_transparent_bg_logo.svg" if Path("assets/openms_transparent_bg_logo.svg").exists() else None,
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


def _get_workspaces_dir() -> Path:
    """Return the workspaces directory from settings.

    Resolved relative to the CWD where ``streamlit run`` was launched, so
    ``"workspaces"`` lands at ``<app-dir>/workspaces/`` for the standard
    invocation. Absolute paths are honored as-is.
    """
    return Path(st.session_state.settings.get("workspaces_dir", "workspaces"))


def _render_workspace_switcher(params: dict) -> None:
    """Render workspace selector, create, and delete controls."""
    workspaces_dir = _get_workspaces_dir()

    with st.expander("**Workspaces**"):
        def change_workspace():
            for key in params:
                if key in st.session_state:
                    del st.session_state[key]
            # Clear page-namespaced keys so stale analysis state from
            # the previous workspace is not shown in the new one.
            stale_keys = [
                k for k in st.session_state
                if k.startswith("peptide_mirror")
            ]
            for key in stale_keys:
                del st.session_state[key]
            st.query_params.workspace = st.session_state["chosen-workspace"]

        workspaces_dir.mkdir(parents=True, exist_ok=True)
        options = sorted(
            f.name for f in workspaces_dir.iterdir() if f.is_dir()
        )
        if not options:
            options = ["default"]

        current_workspace = str(st.session_state.workspace.stem)
        index = options.index(current_workspace) if current_workspace in options else 0

        st.selectbox(
            "choose existing workspace",
            options, index=index,
            on_change=change_workspace, key="chosen-workspace",
        )

        new_ws_name = st.text_input("New workspace name", "", key="new-workspace-name").strip()
        if st.button("**Create Workspace**"):
            if new_ws_name and is_safe_workspace_name(new_ws_name):
                new_path = Path(workspaces_dir, new_ws_name)
                new_path.mkdir(parents=True, exist_ok=True)
                st.session_state.workspace = new_path
                st.query_params.workspace = new_ws_name
                st.rerun()
            else:
                st.warning("Please enter a valid workspace name.")

        if current_workspace != "default":
            confirm_key = "_confirm_delete_workspace"
            if st.session_state.get(confirm_key):
                st.warning(f"Delete workspace **{current_workspace}** and all its data?")
                yes_col, no_col = st.columns(2)
                with yes_col:
                    if st.button("Yes, delete", type="primary", key="ws-confirm-yes"):
                        ws_path = Path(workspaces_dir, current_workspace)
                        if ws_path.exists():
                            shutil.rmtree(ws_path)
                        st.session_state.workspace = Path(workspaces_dir, "default")
                        st.query_params.workspace = "default"
                        del st.session_state[confirm_key]
                        st.rerun()
                with no_col:
                    if st.button("Cancel", key="ws-confirm-no"):
                        del st.session_state[confirm_key]
                        st.rerun()
            else:
                if st.button("Delete this workspace"):
                    st.session_state[confirm_key] = True
                    st.rerun()


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


