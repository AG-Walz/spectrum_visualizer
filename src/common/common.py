"""
Common utilities for MSpecViz Streamlit app.

Adapted from the OpenMS streamlit-template src/common/common.py.
Provides page_setup(), render_sidebar(), workspace management, and UI helpers.
"""

import os
import sys
import json
import time
import shutil

import pandas as pd
import streamlit as st

from typing import Any
from pathlib import Path


# Detect system platform
OS_PLATFORM = sys.platform


def is_safe_workspace_name(name: str) -> bool:
    """
    Check if a workspace name is safe (no path traversal characters).

    Args:
        name: The workspace name to validate.

    Returns:
        bool: True if safe, False if contains path separators or parent references.
    """
    if not name:
        return False
    return "/" not in name and "\\" not in name and name not in ("..", ".")


def load_params(default: bool = False) -> dict[str, Any]:
    """
    Load parameters from a JSON file and return a dictionary containing them.

    If a 'params.json' file exists in the workspace, load the parameters from there.
    Otherwise, load the default parameters from 'settings.json'.

    Args:
        default (bool): Load default parameters. Defaults to False.

    Returns:
        dict[str, Any]: A dictionary containing the parameters.
    """
    if not st.session_state.settings.get("enable_workspaces", False):
        default = True

    path = Path(st.session_state.workspace, "params.json")

    if path.exists() and not default:
        with open(path, "r", encoding="utf-8") as f:
            params = json.load(f)
    else:
        # Use settings.json as default parameters source
        params = dict(st.session_state.settings)

    return params


def save_params(params: dict[str, Any]) -> None:
    """
    Save the given dictionary of parameters to a JSON file in the workspace.

    Captures two layers of state:
    - Top-level settings.json keys already present in params (template pattern).
    - All serializable ``peptide_mirror.*`` session-state keys, stored under a
      ``"peptide_mirror"`` namespace so they can be restored on workspace load.

    Args:
        params (dict[str, Any]): A dictionary containing the parameters to be saved.

    Returns:
        dict[str, Any]: Updated parameters.
    """
    if not st.session_state.settings.get("enable_workspaces", False):
        return

    # Update with any modified session state values (template-level keys)
    for key, value in st.session_state.items():
        if key in params:
            params[key] = value

    # Capture serializable peptide_mirror.* analysis state
    page_state = {}
    for key, value in st.session_state.items():
        if not key.startswith("peptide_mirror."):
            continue
        try:
            json.dumps(value)   # probe serializability; skips VariantGroup dicts, etc.
            page_state[key] = value
        except (TypeError, ValueError):
            pass
    if page_state:
        params["peptide_mirror"] = page_state

    path = Path(st.session_state.workspace, "params.json")
    with open(path, "w", encoding="utf-8") as outfile:
        json.dump(params, outfile, indent=4)

    return params


def page_setup(page: str = "") -> dict[str, Any]:
    """
    Set up the Streamlit page configuration and determine the workspace.

    This function should be run at the start of every page for setup
    and to get the parameters dictionary.

    Args:
        page (str, optional): The name of the current page, by default "".

    Returns:
        dict[str, Any]: A dictionary containing the parameters.
    """
    if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

    # Set Streamlit page configurations
    st.set_page_config(
        page_title=st.session_state.settings.get("app-name", "MSpecViz"),
        page_icon="assets/openms_transparent_bg_logo.svg" if Path("assets/openms_transparent_bg_logo.svg").exists() else None,
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    # Determine the workspace for the current session
    if "workspace" not in st.session_state or (
        "workspace" in st.query_params
        and st.query_params.workspace != st.session_state.workspace.name
    ):
        st.cache_data.clear()
        st.cache_resource.clear()

        st.session_state.location = "local"

        # Define the directory where all workspaces will be stored
        if st.session_state.settings.get("workspaces_dir"):
            workspaces_dir = Path(
                st.session_state.settings["workspaces_dir"],
                "workspaces-" + st.session_state.settings.get("repository-name", "MSpecViz"),
            )
        else:
            workspaces_dir = Path("..")

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

    # Ensure query param is set when workspaces enabled
    if (
        "workspace" not in st.query_params
        and st.session_state.settings.get("enable_workspaces", False)
    ):
        st.query_params.workspace = st.session_state.workspace.name

    # Make sure workspace directories exist
    st.session_state.workspace.mkdir(parents=True, exist_ok=True)

    # Render the sidebar
    params = render_sidebar(page)

    return params


def render_sidebar(page: str = "") -> dict[str, Any]:
    """
    Renders the sidebar on the Streamlit app, which includes the workspace switcher,
    settings, and version info.

    Args:
        page (str): A string indicating the current page.

    Returns:
        dict[str, Any]: The loaded parameters dictionary.
    """
    params = load_params()

    with st.sidebar:
        # Workspace switcher (local mode, workspaces enabled)
        if st.session_state.settings.get("enable_workspaces", False):
            if st.session_state.settings.get("workspaces_dir"):
                workspaces_dir = Path(
                    st.session_state.settings["workspaces_dir"],
                    "workspaces-" + st.session_state.settings.get("repository-name", "MSpecViz"),
                )
            else:
                workspaces_dir = Path("..")

            with st.expander("**Workspaces**"):
                def change_workspace():
                    for key in params:
                        if key in st.session_state:
                            del st.session_state[key]
                    # Clear all page-namespaced keys so stale analysis state from
                    # the previous workspace is not shown in the new one.
                    stale_keys = [
                        k for k in st.session_state
                        if k.startswith("peptide_mirror")
                    ]
                    for key in stale_keys:
                        del st.session_state[key]
                    st.query_params.workspace = st.session_state["chosen-workspace"]

                # Get all available workspaces
                workspaces_dir.mkdir(parents=True, exist_ok=True)
                options = [
                    f.name for f in workspaces_dir.iterdir() if f.is_dir()
                ]
                if not options:
                    options = ["default"]

                current_workspace = str(st.session_state.workspace.stem)
                index = options.index(current_workspace) if current_workspace in options else 0

                st.selectbox(
                    "choose existing workspace",
                    options,
                    index=index,
                    on_change=change_workspace,
                    key="chosen-workspace",
                )

                # Create or remove workspaces
                create_remove = st.text_input("create/remove workspace", "").strip()
                path = Path(workspaces_dir, create_remove)

                if st.button("**Create Workspace**"):
                    if create_remove and is_safe_workspace_name(create_remove):
                        path.mkdir(parents=True, exist_ok=True)
                        st.session_state.workspace = path
                        st.query_params.workspace = create_remove
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.warning("Please enter a valid workspace name.")

                if st.button("Delete Workspace"):
                    if path.exists() and create_remove != "default":
                        shutil.rmtree(path)
                        st.session_state.workspace = Path(workspaces_dir, "default")
                        st.query_params.workspace = "default"
                        st.rerun()

        # Settings expander
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

        # Version display
        with st.container():
            st.markdown(
                """
                <style>
                .version-box {
                    border: 1px solid #a4a5ad;
                    padding: 10px;
                    border-radius: 0.5rem;
                    text-align: center;
                    display: flex;
                    justify-content: center;
                    align-items: center;
                }
                </style>
                """,
                unsafe_allow_html=True,
            )
            version_info = st.session_state.settings.get("version", "0.1.0")
            app_name = st.session_state.settings.get("app-name", "MSpecViz")
            st.markdown(
                f'<div class="version-box">{app_name}<br>Version: {version_info}</div>',
                unsafe_allow_html=True,
            )

    return params


def v_space(n: int, col=None) -> None:
    """
    Writes empty strings to create vertical space in the Streamlit app.

    Args:
        n (int): Number of empty lines to add.
        col: A streamlit column to add vertical space in.
    """
    for _ in range(n):
        if col:
            col.write("")
        else:
            st.write("")


def show_fig(
    fig,
    download_name: str,
    container_width: bool = True,
    selection_session_state_key: str = "",
) -> None:
    """
    Displays a Plotly chart and adds a download button to the plot.

    Args:
        fig: The Plotly figure to display.
        download_name (str): The name for the downloaded file.
        container_width (bool): If True, use container width. Defaults to True.
        selection_session_state_key (str): If set, save selection to session state.
    """
    image_format = st.session_state.get("image-format", "png")

    if not selection_session_state_key:
        st.plotly_chart(
            fig,
            width="stretch" if container_width else "content",
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom", "pan", "select", "lasso",
                    "zoomin", "autoscale", "zoomout", "resetscale",
                ],
                "toImageButtonOptions": {
                    "filename": download_name,
                    "format": image_format,
                },
            },
        )
    else:
        st.plotly_chart(
            fig,
            key=selection_session_state_key,
            selection_mode=["points", "box"],
            on_select="rerun",
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom", "pan", "lasso",
                    "zoomin", "autoscale", "zoomout", "resetscale", "select",
                ],
                "toImageButtonOptions": {
                    "filename": download_name,
                    "format": image_format,
                },
            },
            width="stretch",
        )


def show_table(df: pd.DataFrame, download_name: str = "") -> pd.DataFrame:
    """
    Displays a pandas dataframe and provides a download button.

    Args:
        df (pd.DataFrame): The dataframe to display.
        download_name (str): The name for the downloaded file.

    Returns:
        pd.DataFrame: The displayed dataframe.
    """
    st.dataframe(df, width="stretch")
    if download_name:
        st.download_button(
            "Download Table",
            df.to_csv(sep="\t").encode("utf-8"),
            download_name.replace(" ", "-") + ".tsv",
        )
    return df


def reset_directory(path: Path) -> None:
    """
    Remove the given directory and re-create it.

    Args:
        path (Path): Path to the directory to be reset.
    """
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def get_workspace_url(workspace_name: str) -> str:
    """
    Build a shareable URL for a workspace using the current Streamlit server config.

    Args:
        workspace_name: The workspace name (e.g. "default", "my_run")

    Returns:
        URL string with ?workspace= query parameter
    """
    try:
        import streamlit as st
        # Try to get server address from Streamlit config
        from streamlit.server.server import Server
        port = 8501
    except Exception:
        port = 8501
    return f"?workspace={workspace_name}"


# General warning/error messages
WARNINGS = {
    "missing-mzML": "Upload or select some mzML files first!",
}

ERRORS = {
    "general": "Something went wrong.",
    "workflow": "Something went wrong during workflow execution.",
    "visualization": "Something went wrong during visualization of results.",
}
