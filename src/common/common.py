import json
import os
import shutil
import sys
import uuid
import time
from typing import Any
from pathlib import Path
from streamlit.components.v1 import html

import streamlit as st
import pandas as pd
import psutil
from io import BytesIO
import zipfile

try:
    from tkinter import Tk, filedialog

    TK_AVAILABLE = True
except ImportError:
    TK_AVAILABLE = False

# Detect system platform
OS_PLATFORM = sys.platform


@st.fragment(run_every=5)
def monitor_hardware():
    cpu_progress = psutil.cpu_percent(interval=None) / 100
    ram_progress = 1 - psutil.virtual_memory().available / psutil.virtual_memory().total

    st.text(f"Ram ({ram_progress * 100:.2f}%)")
    st.progress(ram_progress)

    st.text(f"CPU ({cpu_progress * 100:.2f}%)")
    st.progress(cpu_progress)

    st.caption(f"Last fetched at: {time.strftime('%H:%M:%S')}")


def load_params(default: bool = False) -> dict[str, Any]:
    """
    Load parameters from a JSON file and return a dictionary containing them.

    If a 'params.json' file exists in the workspace, load the parameters from there.
    Otherwise, load the default parameters from 'default-parameters.json'.

    Additionally, check if any parameters have been modified by the user during the current session
    and update the values in the parameter dictionary accordingly. Also make sure that all items from
    the parameters dictionary are accessible from the session state as well.

    Args:
        default (bool): Load default parameters. Defaults to True.

    Returns:
        dict[str, Any]: A dictionary containing the parameters.
    """

    # Check if workspace is enabled. If not, load default parameters.
    if not st.session_state.settings.get("enable_workspaces", True):
        default = True

    # Construct the path to the parameter file
    path = Path(st.session_state.workspace, "params.json")

    # Load the parameters from the file, or from the default file if the parameter file does not exist
    if path.exists() and not default:
        with open(path, "r", encoding="utf-8") as f:
            params = json.load(f)
    else:
        # Look for default-parameters.json in the app root directory
        app_root = Path(__file__).parent.parent.parent
        default_params_path = app_root / "default-parameters.json"
        
        if default_params_path.exists():
            with open(default_params_path, "r", encoding="utf-8") as f:
                params = json.load(f)
        else:
            # Fallback to current directory
            with open("default-parameters.json", "r", encoding="utf-8") as f:
                params = json.load(f)

    # Return the parameter dictionary
    return params


def save_params(params: dict[str, Any]) -> None:
    """
    Save the given dictionary of parameters to a JSON file.

    If a 'params.json' file already exists in the workspace, overwrite it with the new parameters.
    Otherwise, create a new 'params.json' file in the workspace directory and save the parameters there.

    Additionally, check if any parameters have been modified by the user during the current session
    and update the values in the parameter dictionary accordingly.

    This function should be run at the end of each page, if the parameters dictionary has been modified directly.
    Note that session states with the same keys will override any direct changes!

    Args:
        params (dict[str, Any]): A dictionary containing the parameters to be saved.

    Returns:
        dict[str, Any]: Updated parameters.
    """

    # Check if the workspace is enabled and if a 'params.json' file exists in the workspace directory
    if not st.session_state.settings.get("enable_workspaces", True):
        return

    # Update the parameter dictionary with any modified parameters from the current session
    for key, value in st.session_state.items():
        if key in params.keys():
            params[key] = value

    # Save the parameter dictionary to a JSON file in the workspace directory
    path = Path(st.session_state.workspace, "params.json")
    with open(path, "w", encoding="utf-8") as outfile:
        json.dump(params, outfile, indent=4)

    return params


def page_setup() -> dict[str, Any]:
    """
    Set up the Streamlit page configuration and determine the workspace for the current session.

    This function should be run at the start of every page for setup and to get the parameters dictionary.

    Returns:
        dict[str, Any]: A dictionary containing the parameters loaded from the parameter file.
    """
    if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

    # Set Streamlit page configurations
    st.set_page_config(
        page_title=st.session_state.settings.get("app-name", "immunuXplore"),
        # page_icon="assets/OpenMS.png",
        # layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )

    # Expand sidebar navigation
    st.markdown(
        """
        <style>
            .stMultiSelect [data-baseweb=select] span{
                max-width: 500px;
                font-size: 1rem;
            }
            div[data-testid='stSidebarNav'] ul {max-height:none}
        </style>
        """,
        unsafe_allow_html=True,
    )

    # st.logo("assets/OpenMS.png")

    # Create google analytics if consent was given
    if (
        ("tracking_consent" not in st.session_state)
        or (st.session_state.tracking_consent is None)
        or (not st.session_state.settings.get("online_deployment", False))
    ):
        st.session_state.tracking_consent = None
    else:
        analytics_settings = st.session_state.settings.get("analytics", {})
        google_analytics = analytics_settings.get("google-analytics", {})
        
        if (google_analytics.get("enabled", False)) and (
            st.session_state.tracking_consent.get("google-analytics", False) == True
        ):
            html(
                """
                <!DOCTYPE html>
                <html lang="en">
                    <head></head>
                    <body><script>
                    window.parent.gtag('consent', 'update', {
                    'analytics_storage': 'granted'
                    });
                    </script></body>
                </html>
                """,
                width=1,
                height=1,
            )
        piwik_settings = analytics_settings.get("piwik-pro", {})
        
        if (piwik_settings.get("enabled", False)) and (
            st.session_state.tracking_consent.get("piwik-pro", False) == True
        ):
            html(
                """
                <!DOCTYPE html>
                <html lang="en">
                    <head></head>
                    <body><script>
                    var consentSettings = {
                        analytics: { status: 1 } // Set Analytics consent to 'on' (1 for on, 0 for off)
                    };
                    window.parent.ppms.cm.api('setComplianceSettings', { consents: consentSettings }, function() {
                        console.log("PiwikPro Analytics consent set to on.");
                    }, function(error) {
                        console.error("Failed to set PiwikPro analytics consent:", error);
                    });
                    </script></body>
                </html>
                """,
                width=1,
                height=1,
            )

    # Determine the workspace for the current session
    if ("workspace" not in st.session_state) or (
        ("workspace" in st.query_params)
        and (st.query_params.workspace != st.session_state.workspace.name)
    ):
        # Clear any previous caches
        st.cache_data.clear()
        st.cache_resource.clear()
        # Check location
        if not st.session_state.settings.get("online_deployment", False):
            st.session_state.location = "local"
            st.session_state["previous_dir"] = os.getcwd()
            st.session_state["local_dir"] = ""
        else:
            st.session_state.location = "online"
        # if we run the packaged windows version, we start within the Python directory -> need to change working directory to ..\st.session_state.settings["repository-name"]
        if "windows" in sys.argv:
            repo_name = st.session_state.settings.get("repository-name", "immunoxplorer-gui")
            os.chdir(Path("..", repo_name))
        # Define the directory where all workspaces will be stored
        workspaces_dir_setting = st.session_state.settings.get("workspaces_dir", "..")
        if (
            workspaces_dir_setting
            and st.session_state.location == "local"
        ):
            repo_name = st.session_state.settings.get("repository-name", "immunoxplorer-gui")
            workspaces_dir = Path(
                workspaces_dir_setting,
                "workspaces-" + repo_name,
            )
        else:
            workspaces_dir = ".."

        # Check if workspace logic is enabled
        if st.session_state.settings.get("enable_workspaces", True):
            if "workspace" in st.query_params:
                st.session_state.workspace = Path(
                    workspaces_dir, st.query_params.workspace
                )
            elif st.session_state.location == "online":
                workspace_id = str(uuid.uuid1())
                st.session_state.workspace = Path(workspaces_dir, workspace_id)
                st.query_params.workspace = workspace_id
            else:
                st.session_state.workspace = Path(workspaces_dir, "default")
                st.query_params.workspace = "default"

        else:
            # Use default workspace when workspace feature is disabled
            st.session_state.workspace = Path(workspaces_dir, "default")


    # If no workspace is specified and workspace feature is enabled, set default workspace and query param
    if (
        "workspace" not in st.query_params
        and st.session_state.settings.get("enable_workspaces", True)
    ):
        st.query_params.workspace = st.session_state.workspace.name

    # Make sure the necessary directories exist
    st.session_state.workspace.mkdir(parents=True, exist_ok=True)

    # Render the sidebar
    params = render_sidebar()

    return params


def render_sidebar() -> None:
    """
    Renders the sidebar on the Streamlit app, which includes the workspace switcher,
    settings, resource utilization & version information.

    Returns:
        None
    """
    params = load_params()
    with st.sidebar:
        # The main page has workspace switcher
        # Display workspace switcher if workspace is enabled in local mode
        if st.session_state.settings.get("enable_workspaces", True):
            with st.expander("🖥️ **Workspaces**"):
                # Workspaces directory specified in the settings.json
                workspaces_dir_setting = st.session_state.settings.get("workspaces_dir", "..")
                if (
                    workspaces_dir_setting
                    and st.session_state.location == "local"
                ):
                    repo_name = st.session_state.settings.get("repository-name", "immunoxplorer-gui")
                    workspaces_dir = Path(
                        workspaces_dir_setting,
                        "workspaces-" + repo_name,
                    )
                else:
                    workspaces_dir = ".."
                # Online: show current workspace name in info text and option to change to other existing workspace
                if st.session_state.location == "local":
                    # Define callback function to change workspace
                    def change_workspace():
                        for key in params.keys():
                            if key in st.session_state.keys():
                                del st.session_state[key]
                        st.session_state.workspace = Path(
                            workspaces_dir, st.session_state["chosen-workspace"]
                        )
                        st.query_params.workspace = st.session_state["chosen-workspace"]

                    # Get all available workspaces as options
                    options = [
                        file.name for file in workspaces_dir.iterdir() if file.is_dir()
                    ]
                    # Let user chose an already existing workspace
                    st.selectbox(
                        "choose existing workspace",
                        options,
                        index=options.index(str(st.session_state.workspace.stem)),
                        on_change=change_workspace,
                        key="chosen-workspace",
                    )
                    # Create or Remove workspaces
                    create_remove = st.text_input("create/remove workspace", "").strip()
                    path = Path(workspaces_dir, create_remove)
                    # Create new workspace
                    if st.button("**Create Workspace**"):
                        if create_remove:
                            path.mkdir(parents=True, exist_ok=True)
                            st.session_state.workspace = path
                            st.query_params.workspace = create_remove
                            # Temporary as the query update takes a short amount of time
                            time.sleep(1)
                            st.rerun()
                        else:
                            st.warning("Please enter a valid workspace name.")
                    # Remove existing workspace and fall back to default
                    if st.button("⚠️ Delete Workspace"):
                        if path.exists():
                            shutil.rmtree(path)
                            st.session_state.workspace = Path(workspaces_dir, "default")
                            st.query_params.workspace = "default"
                            st.rerun()

        # All pages have settings & workflow indicator
        with st.expander("⚙️ **Settings**"):
            img_formats = ["svg", "png", "jpeg", "webp"]
            st.selectbox(
                "image export format",
                img_formats,
                img_formats.index(params["image-format"]),
                key="image-format",
            )

        with st.expander("📊 **Resource Utilization**"):
            monitor_hardware()

        # Display App Version from settings.json
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
            version_info = st.session_state.settings.get("version", "0.0.1")
            app_name = st.session_state.settings.get("app-name", "immunuXplore")
            st.markdown(
                f'<div class="version-box">{app_name}<br>Version: {version_info}</div>',
                unsafe_allow_html=True,
            )
    return params


def v_space(n: int, col=None) -> None:
    """
    Prints empty strings to create vertical space in the Streamlit app.

    Args:
        n (int): An integer representing the number of empty lines to print.
        col: A streamlit column can be passed to add vertical space there.

    Returns:
        None
    """
    for _ in range(n):
        if col:
            col.write("#")
        else:
            st.write("#")


def display_large_dataframe(
    df, chunk_sizes: list[int] = [10, 100, 1_000, 10_000], **kwargs
):
    """
    Displays a large DataFrame in chunks with pagination controls and row selection.

    Args:
        df: The DataFrame to display.
        chunk_sizes: A list of chunk sizes to choose from.
        ...: Additional keyword arguments to pass to the `st.dataframe` function. See: https://docs.streamlit.io/develop/api-reference/data/st.dataframe

    Returns:
        Index of selected row.
    """

    # Dropdown for selecting chunk size
    chunk_size = st.selectbox("Select Number of Rows to Display", chunk_sizes)

    # Calculate total number of chunks
    total_chunks = (len(df) + chunk_size - 1) // chunk_size

    if total_chunks > 1:
        page = int(st.number_input("Select Page", 1, total_chunks, 1, step=1))
    else:
        page = 1

    # Function to get the current chunk of the DataFrame
    def get_current_chunk(df, chunk_size, chunk_index):
        start = chunk_index * chunk_size
        end = min(
            start + chunk_size, len(df)
        )  # Ensure end does not exceed dataframe length
        return df.iloc[start:end], start, end

    # Display the current chunk
    current_chunk_df, start_row, end_row = get_current_chunk(df, chunk_size, page - 1)

    event = st.dataframe(current_chunk_df, **kwargs)

    st.write(
        f"Showing rows {start_row + 1} to {end_row} of {len(df)} ({get_dataframe_mem_useage(current_chunk_df):.2f} MB)"
    )

    rows = event["selection"]["rows"]

    if st.session_state.settings.get("test", False):  # is a test App, return first row as selected
        return 1
    elif not rows:
        return None
    else:
        # Calculate the index based on the current page and chunk size
        base_index = (page - 1) * chunk_size
        print(base_index)
        return base_index + rows[0]


def show_table(df: pd.DataFrame, download_name: str = "") -> None:
    """
    Displays a pandas dataframe using Streamlit's `dataframe` function and
    provides a download button for the same table.

    Args:
        df (pd.DataFrame): The pandas dataframe to display.
        download_name (str): The name to give to the downloaded file. Defaults to empty string.

    Returns:
        df (pd.DataFrame): The possibly edited dataframe.
    """
    # Show dataframe using container width
    st.dataframe(df, use_container_width=True)
    # Show download button with the given download name for the table if name is given
    if download_name:
        st.download_button(
            "Download Table",
            df.to_csv(sep="\t").encode("utf-8"),
            download_name.replace(" ", "-") + ".tsv",
        )
    return df


def show_fig(
    fig,
    download_name: str,
    container_width: bool = True,
    selection_session_state_key: str = "",
) -> None:
    """
    Displays a Plotly chart and adds a download button to the plot.

    Args:
        fig (plotly.graph_objs._figure.Figure): The Plotly figure to display.
        download_name (str): The name for the downloaded file.
        container_width (bool, optional): If True, the figure will use the container width. Defaults to True.
        selection_session_state_key (str, optional): If set, save the rectangular selection to session state with this key.

    Returns:
        None
    """
    if not selection_session_state_key:
        st.plotly_chart(
            fig,
            use_container_width=container_width,
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom",
                    "pan",
                    "select",
                    "lasso",
                    "zoomin",
                    "autoscale",
                    "zoomout",
                    "resetscale",
                ],
                "toImageButtonOptions": {
                    "filename": download_name,
                    "format": st.session_state["image-format"],
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
                    "zoom",
                    "pan",
                    "lasso",
                    "zoomin",
                    "autoscale",
                    "zoomout",
                    "resetscale",
                    "select",
                ],
                "toImageButtonOptions": {
                    "filename": download_name,
                    "format": st.session_state["image-format"],
                },
            },
            use_container_width=True,
        )


def reset_directory(path: Path) -> None:
    """
    Remove the given directory and re-create it.

    Args:
        path (Path): Path to the directory to be reset.

    Returns:
        None
    """
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def get_dataframe_mem_useage(df):
    """
    Get the memory usage of a pandas DataFrame in megabytes.

    Args:
        df (pd.DataFrame): The DataFrame to calculate the memory usage for.

    Returns:
        float: The memory usage of the DataFrame in megabytes.
    """
    # Calculate the memory usage of the DataFrame in bytes
    memory_usage_bytes = df.memory_usage(deep=True).sum()
    # Convert bytes to megabytes
    memory_usage_mb = memory_usage_bytes / (1024**2)
    return memory_usage_mb


def tk_directory_dialog(title: str = "Select Directory", parent_dir: str = os.getcwd()):
    """
    Creates a Tkinter directory dialog for selecting a directory.

    Args:
        title (str): The title of the directory dialog.
        parent_dir (str): The path to the parent directory of the directory dialog.

    Returns:
        str: The path to the selected directory.

    Warning:
        This function is not avaliable in a streamlit cloud context.
    """
    root = Tk()
    root.attributes("-topmost", True)
    root.withdraw()
    file_path = filedialog.askdirectory(title=title, initialdir=parent_dir)
    root.destroy()
    return file_path


def tk_file_dialog(
    title: str = "Select File",
    file_types: list[tuple] = [],
    parent_dir: str = os.getcwd(),
    multiple: bool = True,
):
    """
    Creates a Tkinter file dialog for selecting a file.

    Args:
        title (str): The title of the file dialog.
        file_types (list(tuple)): The file types to filter the file dialog.
        parent_dir (str): The path to the parent directory of the file dialog.
        multiple (bool): If True, multiple files can be selected.

    Returns:
        str: The path to the selected file.

    Warning:
        This function is not avaliable in a streamlit cloud context.
    """
    root = Tk()
    root.attributes("-topmost", True)
    root.withdraw()
    file_types.extend([("All files", "*.*")])
    file_path = filedialog.askopenfilename(
        title=title, filetypes=file_types, initialdir=parent_dir, multiple=True
    )
    root.destroy()
    return file_path


def zip_and_download_files(directory: Path) -> None:
    """
    Creates a zip archive of all files within one or more specified directories,
    including files in subdirectories, and offers it as a download
    button in a Streamlit application.

    Args:
        directory (Path): A directory to be zipped.
    """

    # Flatten all files from all directories (preserving relative paths)
    all_files = []
    for d in directory.iterdir():
        if d.exists() and d.is_dir():
            all_files.extend(list(d.rglob("*")))

    if not all_files:
        st.error("No files to compress.")
        return

    bytes_io = BytesIO()
    n_files = len([f for f in all_files if f.is_file()])
    my_bar = st.progress(0)

    with zipfile.ZipFile(bytes_io, "w", zipfile.ZIP_DEFLATED) as zip_file:
        i = 0
        for file_path in all_files:
            if file_path.is_file():
                # Add file with relative path preserved from parent of base dir
                arcname = file_path.relative_to(directory.parent)
                zip_file.write(file_path, arcname)
                i += 1
                my_bar.progress(i / n_files)

    my_bar.empty()
    bytes_io.seek(0)

    st.download_button(
        label="⬇️ Download Now",
        data=bytes_io,
        file_name=f"{directory.name}.zip",
        mime="application/zip"
    )
    
def load_sample_sheet():
    # loads validated sample sheet from the workspace to session state
    if "sample_sheet" not in st.session_state or st.session_state.sample_sheet is None:
        path = Path(st.session_state["workspace"], "input-valid", "sample-sheet.tsv")
        if path.exists():
            st.session_state.sample_sheet = pd.read_csv(path, sep="\t")
        else:
            st.error("Please upload and validate a dataset first in the **Dataset Selection** section.")
            st.stop()