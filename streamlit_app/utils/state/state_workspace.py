"""Workspace directory helpers (name validation + path resolution) used by state_sidebar and page_setup."""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from streamlit_app.utils.state.state_core import PROJECT_ROOT


def is_safe_workspace_name(name: str) -> bool:
    """Reject empty names and any name containing path-traversal characters (/, \\, ., ..)."""
    if not name:
        return False
    return "/" not in name and "\\" not in name and name not in ("..", ".")


def _get_workspaces_dir() -> Path:
    """Return the workspaces root from settings.workspaces_dir — relative paths are anchored to the project root, not CWD."""
    configured = Path(st.session_state.settings.get("workspaces_dir", "workspaces"))
    if configured.is_absolute():
        return configured
    return PROJECT_ROOT / configured
