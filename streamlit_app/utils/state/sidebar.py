"""
Sidebar, workspace switcher, and page-setup entry point.

``page_setup()`` is called at the top of every view; ``render_sidebar``
renders the shared left-hand controls; ``_render_workspace_switcher`` is the
workspace selector within the sidebar.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import streamlit as st

from mspecviz.common.defaults import DEFAULTS

from streamlit_app.utils.state.core import PROJECT_ROOT, STATE_PREFIX, _persisted_key
from streamlit_app.utils.state.params import load_params, save_params
from streamlit_app.utils.state.workspace import (
    _get_workspaces_dir,
    is_safe_workspace_name,
)


def _render_workspace_switcher(params: dict) -> None:
    """Render workspace selector, create, and delete controls."""
    workspaces_dir = _get_workspaces_dir()

    with st.expander("**Workspaces**"):
        def change_workspace():
            # Flush any in-flight session state to the OUTGOING workspace's
            # params.json BEFORE deleting keys. Without this, peptides the
            # user typed in a text_area whose blur event races with the
            # dropdown on_change can be lost on switch.
            save_params(params)

            for key in params:
                st.session_state.pop(key, None)
            # Clear page-namespaced keys so stale analysis state from
            # the previous workspace is not shown in the new one.
            stale_keys = [
                k for k in st.session_state
                if k.startswith(STATE_PREFIX)
            ]
            for key in stale_keys:
                st.session_state.pop(key, None)
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
                        st.session_state.pop(confirm_key, None)
                        st.rerun()
                with no_col:
                    if st.button("Cancel", key="ws-confirm-no"):
                        st.session_state.pop(confirm_key, None)
                        st.rerun()
            else:
                if st.button("Delete this workspace"):
                    st.session_state[confirm_key] = True
                    st.rerun()


def page_setup(page: str = "") -> dict[str, Any]:
    """Set up page config, workspace, and sidebar. Call at the top of every page."""
    if "settings" not in st.session_state:
        settings_path = PROJECT_ROOT / "streamlit_app" / "config" / "settings.json"
        with open(settings_path, "r") as f:
            st.session_state.settings = json.load(f)

    logo_path = PROJECT_ROOT / "assets" / "openms_transparent_bg_logo.svg"
    st.set_page_config(
        page_title=st.session_state.settings.get("app_name", "MSpecViz"),
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
            persisted = _persisted_key("image_format")
            # Seed session_state BEFORE the widget is created so Streamlit
            # uses the restored value on page navigation. Reading from the
            # non-widget persisted backup survives widget garbage collection.
            if "image_format" not in st.session_state:
                restored = st.session_state.get(
                    persisted, params.get("image_format", DEFAULTS.image_format)
                )
                if restored not in img_formats:
                    restored = DEFAULTS.image_format
                st.session_state["image_format"] = restored
            selected_format = st.selectbox(
                "image export format",
                img_formats,
                key="image_format",
            )
            st.session_state[persisted] = selected_format
            params["image_format"] = selected_format
            save_params(params)

        version_info = st.session_state.settings.get("version", "0.1.0")
        app_name = st.session_state.settings.get("app_name", "MSpecViz")
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
