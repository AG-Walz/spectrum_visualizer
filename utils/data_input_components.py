"""Shared data input UI components for mzML and idXML file selection.

This module provides the unified file selection interface used by both:
- spectrum_analysis.py: Select files for spectral comparison
- dataset_filter.py: Select files for filtering operations

DESIGN:
Two independent steps:
1. Server Selection: Choose data source (Local Computer / Current Server / SSH servers)
2. File Selection: Choose method (Individual files / Folder browsing)

Both file selection methods work for ALL server types:
- Local: Upload files individually OR browse local folders
- Server: Enter paths individually OR browse server folders

IMPORTANT: This is the single source of truth for file input UI.
Edit once here, affects all workspaces.
"""

import streamlit as st
import os
import socket
import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# Import shared utility functions
from file_selection_components import (
    render_directory_browser,
    get_files_with_substring
)

from shared_imports import (
    create_file_options_with_sizes,
    map_display_names_to_paths
)

logger = logging.getLogger(__name__)


def convert_uploaded_file_to_path(file_obj) -> str:
    """Convert UploadedFile object to temporary file path.
    
    For UploadedFile objects, saves to temp directory only when needed.
    For path strings, returns as-is.
    
    Args:
        file_obj: Either an UploadedFile object or a file path string
        
    Returns:
        File path string
    """
    if hasattr(file_obj, 'read') and hasattr(file_obj, 'name'):
        # It's an UploadedFile - save to temp directory
        temp_dir = tempfile.gettempdir()
        temp_path = os.path.join(temp_dir, file_obj.name)
        
        # Save the uploaded file
        with open(temp_path, 'wb') as f:
            f.write(file_obj.getvalue())
        
        logger.info(f"Saved uploaded file to: {temp_path}")
        return temp_path
    else:
        # Already a path string
        return str(file_obj)


def get_available_servers() -> Dict[str, Dict[str, Any]]:
    """
    Get list of available servers for file access.
    Returns dict of server configurations with display names.
    """
    servers = {}
    
    # Local computer option (for file upload)
    servers["Local Computer"] = {
        "type": "local",
        "path": None,
        "display": "💻 Local Computer",
        "description": "Upload files from your local machine"
    }
    
    # Try to read SSH config for remote servers
    ssh_config_path = Path.home() / ".ssh" / "config"
    if ssh_config_path.exists():
        try:
            with open(ssh_config_path, 'r') as f:
                config_lines = f.readlines()
                
            current_host = None
            current_config = {}
            
            for line in config_lines:
                line = line.strip()
                if line.startswith("Host ") and not line.startswith("Host *"):
                    # Save previous host config
                    if current_host and current_host not in ["*", "github.com"]:
                        servers[current_host] = {
                            "type": "ssh",
                            "path": current_config.get("path", "/mnt/volume/"),
                            "hostname": current_config.get("hostname", current_host),
                            "user": current_config.get("user", os.getenv("USER", "user")),
                            "display": f"🖥️ {current_host}",
                            "description": f"SSH: {current_config.get('hostname', current_host)}"
                        }
                    
                    # Start new host
                    current_host = line.split()[1]
                    current_config = {}
                
                elif current_host:
                    if line.startswith("HostName "):
                        current_config["hostname"] = line.split()[1]
                    elif line.startswith("User "):
                        current_config["user"] = line.split()[1]
            
            # Save last host
            if current_host and current_host not in ["*", "github.com"]:
                servers[current_host] = {
                    "type": "ssh",
                    "path": current_config.get("path", "/mnt/volume/"),
                    "hostname": current_config.get("hostname", current_host),
                    "user": current_config.get("user", os.getenv("USER", "user")),
                    "display": f"🖥️ {current_host}",
                    "description": f"SSH: {current_config.get('hostname', current_host)}"
                }
        
        except Exception as e:
            logger.warning(f"Could not read SSH config from {ssh_config_path}: {e}")
    
    # Add current server (where script is running)
    current_hostname = socket.gethostname()
    if "Current Server" not in servers:
        servers["Current Server"] = {
            "type": "current",
            "path": "/mnt/volume/",  # Default from settings
            "display": f"🌐 Current Server ({current_hostname})",
            "description": f"Server where this script is running: {current_hostname}"
        }
    
    return servers


def render_server_directory_selector(default_path: Optional[str] = None, session_key: str = 'server_default', server_name: Optional[str] = None):
    """Render a directory selector for server files with folder browsing.
    
    Args:
        default_path: Default directory path to start from
        session_key: Base session key for the directory browser
        server_name: Server name to append to session key (for server-specific state)
    """
    if default_path is None:
        default_path = st.session_state.settings.get("data_directories", {}).get("default", "/mnt/volume/")
    
    # Use provided server_name or fall back to session state
    if server_name is None:
        server_name = st.session_state.get(f'{session_key}_selected_server_name', 'default')
    
    return render_directory_browser(
        session_key=f'{session_key}_dir_{server_name}',
        default_path=default_path,
        show_subdirs=True,
        compact=True
    )


def render_unified_file_browser(source: str = "server", server_info: Optional[Dict] = None, session_key_prefix: str = "unified"):
    """
    Unified file browser interface that works for all server types.
    Shows both file selection methods (individual + folder browsing) regardless of source.
    Files from both tabs are merged and accumulated in session state.
    
    Args:
        source: "server" for remote/current server files, "local" for local uploads
        server_info: Dictionary with server connection information
        session_key_prefix: Prefix for session state keys to avoid conflicts between workspaces
    
    Returns:
        (mzml_files, idxml_files) - Lists of file paths (UploadedFile objects converted to temp paths)
    """
    # Initialize session state for accumulated files
    mzml_session_key = f"{session_key_prefix}_accumulated_mzml"
    idxml_session_key = f"{session_key_prefix}_accumulated_idxml"
    
    if mzml_session_key not in st.session_state:
        st.session_state[mzml_session_key] = []
    if idxml_session_key not in st.session_state:
        st.session_state[idxml_session_key] = []
    
    # Collect files from current selection
    mzml_files = []
    idxml_files = []
    
    # Add file upload option for local computer
    if source == "local":
        st.markdown("### 📤 Option 1: Upload Files from Browser")
        st.caption("Upload files directly from wherever your browser is running")
        
        col_upload1, col_upload2 = st.columns(2)
        
        with col_upload1:
            st.markdown("**Upload mzML Files**")
            uploaded_mzml = st.file_uploader(
                "Choose mzML files",
                type=['mzML', 'mzml'],
                accept_multiple_files=True,
                key=f"{session_key_prefix}_upload_mzml",
                help="Upload mzML files from your local machine"
            )
            if uploaded_mzml:
                total_mb = sum(f.size for f in uploaded_mzml) / 1024 / 1024
                st.info(f"📊 {len(uploaded_mzml)} file(s) selected ({total_mb:.1f} MB)")
                if st.button("➕ Add Uploaded mzML Files", key=f"{session_key_prefix}_add_uploaded_mzml", type="primary"):
                    with st.spinner("Saving mzML files..."):
                        for uploaded_file in uploaded_mzml:
                            # Convert to temp path and add to session state
                            temp_path = convert_uploaded_file_to_path(uploaded_file)
                            if temp_path not in st.session_state[mzml_session_key]:
                                st.session_state[mzml_session_key].append(temp_path)
                    st.success(f"✅ Added {len(uploaded_mzml)} mzML file(s)")
                    st.rerun()
        
        with col_upload2:
            st.markdown("**Upload idXML Files**")
            uploaded_idxml = st.file_uploader(
                "Choose idXML files",
                type=['idXML', 'idxml', 'mzid', 'mzidentml'],
                accept_multiple_files=True,
                key=f"{session_key_prefix}_upload_idxml",
                help="Upload idXML/mzIdentML files from your local machine"
            )
            if uploaded_idxml:
                total_mb = sum(f.size for f in uploaded_idxml) / 1024 / 1024
                st.info(f"📊 {len(uploaded_idxml)} file(s) selected ({total_mb:.1f} MB)")
                if st.button("➕ Add Uploaded idXML Files", key=f"{session_key_prefix}_add_uploaded_idxml", type="primary"):
                    with st.spinner("Saving idXML files..."):
                        for uploaded_file in uploaded_idxml:
                            # Convert to temp path and add to session state
                            temp_path = convert_uploaded_file_to_path(uploaded_file)
                            if temp_path not in st.session_state[idxml_session_key]:
                                st.session_state[idxml_session_key].append(temp_path)
                    st.success(f"✅ Added {len(uploaded_idxml)} idXML file(s)")
                    st.rerun()
        
        st.markdown("---")
        st.markdown("### 📁 Option 2: Browse Local Folders")
        st.caption("Browse directories accessible from the server")
    else:
        st.markdown("### 📁 Browse Server Folders")
        st.caption("💡 You can add multiple folders - all files from all folders will be combined for analysis")
    
    # Set default path based on source
    if source == "local":
        # For local computer, start from user's home directory
        default_path = str(Path.home())
    else:
        # For servers, use the configured server path
        default_path = server_info.get("path", "/mnt/volume/") if server_info else "/mnt/volume/"
    
    # Allow multiple folder selection (use session_key_prefix for isolation)
    folders_key = f'{session_key_prefix}_selected_folders'
    if folders_key not in st.session_state:
        st.session_state[folders_key] = []
    
    col_a, col_b = st.columns([4, 1])
    with col_a:
        # Get current server name from session state
        current_server_name = st.session_state.get(f'{session_key_prefix}_selected_server_name', 'default')
        current_folder = render_server_directory_selector(
            default_path=default_path, 
            session_key=session_key_prefix,
            server_name=current_server_name
        )
    with col_b:
        # Add spacing to align with directory browser (accounts for back button row)
        st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
        if current_folder and st.button("➕ Add Folder", key=f"{session_key_prefix}_add_folder_btn", use_container_width=True):
            if current_folder not in st.session_state[folders_key]:
                st.session_state[folders_key].append(current_folder)
                st.rerun()
    
    # Display selected folders with remove option in a collapsible section
    if st.session_state[folders_key]:
        with st.expander(f"📁 Selected Folders ({len(st.session_state[folders_key])})", expanded=False):
            folders_to_remove = []
            for idx, folder in enumerate(st.session_state[folders_key]):
                col_folder, col_remove = st.columns([9, 1])
                with col_folder:
                    st.text(folder)
                with col_remove:
                    if st.button("✕", key=f"{session_key_prefix}_remove_folder_{idx}", help="Remove this folder"):
                        folders_to_remove.append(folder)
            
            # Remove folders after iteration
            for folder in folders_to_remove:
                st.session_state[folders_key].remove(folder)
                st.rerun()
    
    selected_folders = st.session_state[folders_key]
    
    if selected_folders:
    
        # Show filtering options
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**📊 mzML Files**")
            mzml_substring = st.text_input(
                "Filter mzML files (optional)",
                value="",
                key=f"{session_key_prefix}_mzml_filter",
                placeholder="e.g., HNSCC, Rep1",
                help="Filter by substring in filename"
            )
            
            mzml_include_subdirs = st.checkbox(
                "Include subdirectories",
                value=True,
                key=f"{session_key_prefix}_mzml_subdirs"
            )
        
        with col2:
            st.markdown("**🔍 idXML Files**")
            idxml_substring = st.text_input(
                "Filter idXML files (optional)",
                value="",
                key=f"{session_key_prefix}_idxml_filter",
                placeholder="e.g., HNSCC, Rep1",
                help="Filter by substring in filename"
            )
            
            idxml_include_subdirs = st.checkbox(
                "Include subdirectories",
                value=True,
                key=f"{session_key_prefix}_idxml_subdirs"
            )
        
        # Find files with filters from all selected folders
        mzml_files_found = []
        idxml_files_found = []
        
        for folder in selected_folders:
            folder_mzml = get_files_with_substring(
                folder, 
                mzml_substring, 
                ['.mzML']
            )
            folder_idxml = get_files_with_substring(
                folder,
                idxml_substring,
                ['.idXML', '.mzid', '.mzidentml']
            )
            mzml_files_found.extend(folder_mzml)
            idxml_files_found.extend(folder_idxml)
        
        # Remove duplicates while preserving order
        mzml_files_found = list(dict.fromkeys(mzml_files_found))
        idxml_files_found = list(dict.fromkeys(idxml_files_found))
        
        # Show results with folder breakdown
        if len(selected_folders) > 1:
            with st.expander(f"📊 Files by Folder (Total: {len(mzml_files_found)} mzML, {len(idxml_files_found)} idXML)", expanded=False):
                for folder in selected_folders:
                    folder_mzml_count = len([f for f in mzml_files_found if f.startswith(folder)])
                    folder_idxml_count = len([f for f in idxml_files_found if f.startswith(folder)])
                    st.text(f"📁 {folder}")
                    st.text(f"   └─ {folder_mzml_count} mzML, {folder_idxml_count} idXML")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if mzml_files_found:
                count_msg = f"Found {len(mzml_files_found)} mzML files across {len(selected_folders)} folder(s)"
                if len(mzml_files_found) >= 10000:
                    st.warning(f"⚠️ {count_msg} (search limit reached - results may be incomplete)")
                else:
                    st.info(count_msg)
                
                select_all_mzml = st.checkbox(
                    f"Select all {len(mzml_files_found)} mzML files",
                    key=f"{session_key_prefix}_select_all_mzml"
                )
                
                if select_all_mzml:
                    mzml_files = mzml_files_found
                else:
                    file_options = create_file_options_with_sizes(mzml_files_found)
                    selected_display = st.multiselect(
                        "Select mzML files",
                        options=file_options,
                        key=f"{session_key_prefix}_mzml_select"
                    )
                    mzml_files = map_display_names_to_paths(selected_display, mzml_files_found)
            else:
                st.warning("No mzML files found")
        
        with col2:
            if idxml_files_found:
                count_msg = f"Found {len(idxml_files_found)} idXML files across {len(selected_folders)} folder(s)"
                if len(idxml_files_found) >= 10000:
                    st.warning(f"⚠️ {count_msg} (search limit reached - results may be incomplete)")
                else:
                    st.info(count_msg)
                
                select_all_idxml = st.checkbox(
                    f"Select all {len(idxml_files_found)} idXML files",
                    key=f"{session_key_prefix}_select_all_idxml"
                )
                
                if select_all_idxml:
                    idxml_files = idxml_files_found
                else:
                    file_options = create_file_options_with_sizes(idxml_files_found)
                    selected_display = st.multiselect(
                        "Select idXML files",
                        options=file_options,
                        key=f"{session_key_prefix}_idxml_select"
                    )
                    idxml_files = map_display_names_to_paths(selected_display, idxml_files_found)
            else:
                st.warning("No idXML files found")
        
        # Add button to accumulate folder selections
        if mzml_files or idxml_files:
            if st.button("➕ Add These Files to Selection", key=f"{session_key_prefix}_add_folder_files", type="primary"):
                for f in mzml_files:
                    if f not in st.session_state[mzml_session_key]:
                        st.session_state[mzml_session_key].append(f)
                for f in idxml_files:
                    if f not in st.session_state[idxml_session_key]:
                        st.session_state[idxml_session_key].append(f)
                st.success(f"✅ Added {len(mzml_files)} mzML + {len(idxml_files)} idXML files")
            st.rerun()
    
    # Return accumulated files from session state
    return st.session_state[mzml_session_key], st.session_state[idxml_session_key]


def render_data_input_section(session_key_prefix: str = "data_input", show_process_button: bool = True, process_callback = None):
    """
    Render the unified data input section with independent server and file selection.
    
    This is the main entry point for file selection UI used by all workspaces.
    Server selection is done first and saved independently. File selection method
    (individual vs folder browsing) is then available regardless of server choice.
    
    Args:
        session_key_prefix: Unique prefix for session state keys (prevents conflicts between workspaces)
        show_process_button: Whether to show the "Confirm File Selection" button
        process_callback: Optional callback function(mzml_files, idxml_files) called when files are confirmed
    
    Returns:
        (mzml_files, idxml_files) - Lists of selected files
    """
    st.header("📁 Data Input")
    
    # Get available servers
    with st.spinner("Loading server list..."):
        available_servers = get_available_servers()
    
    # ========================================================================
    # STEP 1: Server Selection (Independent)
    # ========================================================================
    st.subheader("🌐 Step 1: Select Server")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        # Create server options with display names
        server_options = [info["display"] for info in available_servers.values()]
        server_names = list(available_servers.keys())
        
        selected_display = st.selectbox(
            "Choose data source:",
            options=server_options,
            key=f"{session_key_prefix}_server_selector",
            help="Select where your files are located - this choice is independent of how you select files below"
        )
        
        # Map back to server name
        selected_server_name = server_names[server_options.index(selected_display)]
        selected_server = available_servers[selected_server_name]
        
        # Show server description
        st.caption(selected_server["description"])
    
    with col2:
        if st.button("🔄 Reset", key=f"{session_key_prefix}_reset", help="Clear all selections and start over"):
            # Clear relevant session state for this workspace
            keys_to_clear = [
                f'{session_key_prefix}_selected_folders',
                f'{session_key_prefix}_selected_server',
                f'{session_key_prefix}_selected_server_name',
                f'{session_key_prefix}_individual_mzml_upload',
                f'{session_key_prefix}_individual_idxml_upload',
                f'{session_key_prefix}_individual_mzml_paths',
                f'{session_key_prefix}_individual_idxml_paths',
            ]
            for key in keys_to_clear:
                if key in st.session_state:
                    del st.session_state[key]
            st.rerun()
    
    # Check if server changed and clear directory state if so
    previous_server = st.session_state.get(f'{session_key_prefix}_selected_server_name')
    if previous_server and previous_server != selected_server_name:
        # Clear directory browser state when switching servers
        old_dir_key = f'{session_key_prefix}_dir_{previous_server}'
        if old_dir_key in st.session_state:
            del st.session_state[old_dir_key]
        # Also clear selected folders when switching servers
        folders_key = f'{session_key_prefix}_selected_folders'
        if folders_key in st.session_state:
            st.session_state[folders_key] = []
    
    # Store selected server in session state (with prefix for isolation)
    st.session_state[f'{session_key_prefix}_selected_server'] = selected_server
    st.session_state[f'{session_key_prefix}_selected_server_name'] = selected_server_name
    
    st.markdown("---")
    
    # ========================================================================
    # STEP 2: File Selection (Independent of server choice)
    # ========================================================================
    st.subheader("📂 Step 2: Select Files")
    
    server_type = selected_server["type"]
    source = "local" if server_type == "local" else "server"
    
    # Show accumulated files summary
    mzml_session_key = f"{session_key_prefix}_accumulated_mzml"
    idxml_session_key = f"{session_key_prefix}_accumulated_idxml"
    
    if mzml_session_key in st.session_state or idxml_session_key in st.session_state:
        mzml_count = len(st.session_state.get(mzml_session_key, []))
        idxml_count = len(st.session_state.get(idxml_session_key, []))
        
        if mzml_count > 0 or idxml_count > 0:
            col_summary, col_clear = st.columns([3, 1])
            with col_summary:
                st.success(f"📋 **Selected:** {mzml_count} mzML files, {idxml_count} idXML files")
            with col_clear:
                if st.button("🗑️ Clear All", key=f"{session_key_prefix}_clear_all", help="Remove all selected files"):
                    st.session_state[mzml_session_key] = []
                    st.session_state[idxml_session_key] = []
                    st.rerun()
            
            # Show file lists in expander
            with st.expander("📂 View Selected Files", expanded=False):
                if mzml_count > 0:
                    st.markdown(f"**mzML Files ({mzml_count}):**")
                    for i, f in enumerate(st.session_state[mzml_session_key][:10], 1):
                        st.text(f"{i}. {Path(f).name}")
                    if mzml_count > 10:
                        st.text(f"... and {mzml_count - 10} more")
                
                if idxml_count > 0:
                    st.markdown(f"**idXML Files ({idxml_count}):**")
                    for i, f in enumerate(st.session_state[idxml_session_key][:10], 1):
                        st.text(f"{i}. {Path(f).name}")
                    if idxml_count > 10:
                        st.text(f"... and {idxml_count - 10} more")
    
    st.markdown("---")
    
    # Render unified file browser (same two-tab interface for all sources)
    mzml_files, idxml_files = render_unified_file_browser(
        source=source, 
        server_info=selected_server,
        session_key_prefix=session_key_prefix
    )
    
    # Process files if callback provided and button shown
    if show_process_button and mzml_files and idxml_files:
        if process_callback:
            process_callback(mzml_files, idxml_files)
    
    return mzml_files, idxml_files