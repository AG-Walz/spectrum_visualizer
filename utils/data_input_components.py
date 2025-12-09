"""Shared data input UI components for mzML and idXML file selection.

This module provides consistent file selection interfaces used by:
- spectrum_analysis.py: Select files for spectral comparison
- dataset_filter.py: Select files for filtering operations

The data input flow is separated into two independent steps:
1. Server/Source Selection (where files are located)
2. File Selection (which specific files to use)

Usage:
    from data_input_components import render_complete_data_input_section
    
    # Complete workflow (server + files)
    mzml_files, idxml_files, base_dir = render_complete_data_input_section(...)
    
    # Or use individual components
    selected_server = render_server_selection_section()
    mzml_files, idxml_files = render_file_selection_section(server_info)
"""

import streamlit as st
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
import socket

# Import shared components
from file_selection_components import (
    render_directory_browser,
    render_file_filter_and_selection
)


def render_mzml_idxml_file_selector(
    session_key_prefix: str = "data_input",
    require_both: bool = False,
    show_mzml: bool = True,
    show_idxml: bool = True,
    header_text: Optional[str] = None,
    description: Optional[str] = None
) -> Tuple[List[str], List[str], str]:
    """
    Render unified file selection interface for mzML and idXML files.
    
    This provides a consistent UI across all workspaces:
    1. Directory browser to navigate filesystem
    2. Side-by-side mzML and idXML file selection with filtering
    3. Selection summary
    
    Args:
        session_key_prefix: Prefix for session state keys (ensures isolation between workspaces)
        require_both: If True, shows warning when only one file type is selected
        show_mzml: Whether to show mzML file selection column
        show_idxml: Whether to show idXML file selection column
        header_text: Custom header text (default: "📁 Select Data Files")
        description: Optional description text shown below header
    
    Returns:
        (selected_mzml_paths, selected_idxml_paths, current_directory)
    """
    
    # Header
    if header_text:
        st.header(header_text)
    else:
        st.header("📁 Select Data Files")
    
    if description:
        st.markdown(description)
    
    # Get default directory from settings
    default_dir = st.session_state.settings.get("data_directories", {}).get("default", "/mnt/volume/")
    
    # Directory browser (shared across both file types)
    current_dir = render_directory_browser(
        session_key=f'{session_key_prefix}_browser',
        default_path=default_dir,
        show_subdirs=True,
        compact=True
    )
    
    st.markdown("---")
    
    # Determine column layout
    if show_mzml and show_idxml:
        col1, col2 = st.columns(2)
    elif show_mzml:
        col1 = st.container()
        col2 = None
    elif show_idxml:
        col1 = None
        col2 = st.container()
    else:
        st.error("⚠️ Must show at least one file type (mzML or idXML)")
        return [], [], current_dir
    
    selected_mzml = []
    selected_idxml = []
    
    # mzML file selection
    if show_mzml and col1:
        with col1:
            st.subheader("📊 mzML Files")
            st.caption("Mass spectrometry data files")
            
            selected_mzml, _, _ = render_file_filter_and_selection(
                directory=current_dir,
                file_type="mzML",
                extensions=['.mzml', '.mzML'],
                session_key=f"{session_key_prefix}_mzml",
                default_substring="",
                include_subdirs=True,
                show_select_all=True,
                multiselect=True
            )
    
    # idXML file selection
    if show_idxml and col2:
        with col2:
            st.subheader("🔍 idXML Files")
            st.caption("Peptide identification files")
            
            selected_idxml, _, _ = render_file_filter_and_selection(
                directory=current_dir,
                file_type="idXML",
                extensions=['.idxml', '.idXML', '.mzid', '.mzidentml'],
                session_key=f"{session_key_prefix}_idxml",
                default_substring="",
                include_subdirs=True,
                show_select_all=True,
                multiselect=True
            )
    
    # Selection summary
    st.markdown("---")
    
    if selected_mzml and selected_idxml:
        st.success(f"✅ **Selected:** {len(selected_mzml)} mzML file(s) and {len(selected_idxml)} idXML file(s)")
    elif selected_mzml:
        if require_both:
            st.warning(f"⚠️ Selected {len(selected_mzml)} mzML file(s) — idXML files also required")
        else:
            st.info(f"📊 Selected {len(selected_mzml)} mzML file(s)")
    elif selected_idxml:
        if require_both:
            st.warning(f"⚠️ Selected {len(selected_idxml)} idXML file(s) — mzML files also required")
        else:
            st.info(f"🔍 Selected {len(selected_idxml)} idXML file(s)")
    else:
        st.info("👆 Select files from the directory above to continue")
    
    return selected_mzml, selected_idxml, current_dir


def render_idxml_only_selector(
    session_key_prefix: str = "idxml_input",
    header_text: Optional[str] = None,
    description: Optional[str] = None
) -> Tuple[List[str], str]:
    """
    Simplified interface for idXML-only selection (wrapper for consistency).
    
    Args:
        session_key_prefix: Prefix for session state keys
        header_text: Custom header text
        description: Optional description text
    
    Returns:
        (selected_idxml_paths, current_directory)
    """
    _, selected_idxml, current_dir = render_mzml_idxml_file_selector(
        session_key_prefix=session_key_prefix,
        require_both=False,
        show_mzml=False,
        show_idxml=True,
        header_text=header_text or "📁 Select idXML Files",
        description=description
    )
    
    return selected_idxml, current_dir


def render_mzml_only_selector(
    session_key_prefix: str = "mzml_input",
    header_text: Optional[str] = None,
    description: Optional[str] = None
) -> Tuple[List[str], str]:
    """
    Simplified interface for mzML-only selection (wrapper for consistency).
    
    Args:
        session_key_prefix: Prefix for session state keys
        header_text: Custom header text
        description: Optional description text
    
    Returns:
        (selected_mzml_paths, current_directory)
    """
    selected_mzml, _, current_dir = render_mzml_idxml_file_selector(
        session_key_prefix=session_key_prefix,
        require_both=False,
        show_mzml=True,
        show_idxml=False,
        header_text=header_text or "📁 Select mzML Files",
        description=description
    )
    
    return selected_mzml, current_dir


# ================================================================
# COMPLETE DATA INPUT WORKFLOW (SERVER + FILE SELECTION)
# ================================================================

def get_available_servers() -> Dict[str, Dict[str, Any]]:
    """
    Get list of available servers/data sources.
    
    Returns:
        Dictionary of server configurations:
        {
            'server_key': {
                'display': 'Display Name',
                'description': 'Server description',
                'type': 'current'|'local'|'ssh',
                'path': '/default/path'
            }
        }
    """
    servers = {}
    
    # Current server (always available)
    hostname = socket.gethostname()
    servers['current'] = {
        'display': f'📡 Current Server ({hostname})',
        'description': f'Files on this server: {hostname}',
        'type': 'current',
        'path': '/mnt/volume/'
    }
    
    # Local computer upload
    servers['local'] = {
        'display': '💻 Local Computer',
        'description': 'Upload files from your local computer',
        'type': 'local',
        'path': None
    }
    
    # TODO: Add SSH servers from config if available
    # This could read from ~/.ssh/config or a custom configuration file
    
    return servers


def render_server_selection_section(
    session_key: str = "server_selection",
    header_text: Optional[str] = None
) -> Dict[str, Any]:
    """
    Render server/source selection interface (Step 1).
    
    This is the first step in the data input workflow, separate from file selection.
    User chooses WHERE the data is located (current server, local computer, remote server).
    
    Args:
        session_key: Session state key for storing selection
        header_text: Custom header text
    
    Returns:
        Selected server information dictionary
    """
    if header_text:
        st.header(header_text)
    else:
        st.header("🌐 Step 1: Select Data Source")
    
    st.caption("Choose where your data files are located")
    
    # Get available servers
    available_servers = get_available_servers()
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        # Create server options with display names
        server_options = [info["display"] for info in available_servers.values()]
        server_names = list(available_servers.keys())
        
        # Check if we have a previous selection
        previous_selection = st.session_state.get(f'{session_key}_name', 'current')
        default_index = server_names.index(previous_selection) if previous_selection in server_names else 0
        
        selected_display = st.selectbox(
            "Data source:",
            options=server_options,
            index=default_index,
            key=f"{session_key}_selector",
            help="Select where your mzML and idXML files are located"
        )
        
        # Map back to server name
        selected_server_name = server_names[server_options.index(selected_display)]
        selected_server = available_servers[selected_server_name]
        
        # Show server description
        st.caption(f"ℹ️ {selected_server['description']}")
    
    with col2:
        st.write("")  # Spacing
        if st.button("🔄 Reset All", key=f"{session_key}_reset", help="Clear all selections and start over"):
            # Clear all data input related session state
            keys_to_clear = [k for k in st.session_state.keys() if 
                           any(x in k for x in ['_mzml', '_idxml', '_browser', '_filter', 
                                               'selected_', 'current_dir', session_key])]
            for key in keys_to_clear:
                del st.session_state[key]
            st.rerun()
    
    # Store selected server in session state
    st.session_state[f'{session_key}_name'] = selected_server_name
    st.session_state[f'{session_key}_info'] = selected_server
    
    # Show warning for local uploads if needed
    if selected_server['type'] == 'local':
        st.info("💡 **Local Upload Mode**: You'll upload files from your computer in the next step")
    elif selected_server['type'] == 'ssh':
        st.warning("⚠️ **SSH Server**: Requires files to be mounted or accessible via filesystem")
    
    return selected_server


def render_file_selection_section(
    server_info: Dict[str, Any],
    session_key_prefix: str = "file_selection",
    require_both: bool = False,
    show_mzml: bool = True,
    show_idxml: bool = True,
    header_text: Optional[str] = None,
    allow_upload: bool = True
) -> Tuple[List[str], List[str], str]:
    """
    Render file selection interface (Step 2) based on selected server.
    
    This is the second step, independent from server selection.
    Shows appropriate file selection UI based on the server type:
    - Local: File upload interface
    - Current/SSH: File browser interface
    
    Args:
        server_info: Server configuration from render_server_selection_section()
        session_key_prefix: Prefix for session state keys
        require_both: Whether both mzML and idXML are required
        show_mzml: Show mzML selection
        show_idxml: Show idXML selection
        header_text: Custom header text
        allow_upload: For local mode, whether to allow file uploads
    
    Returns:
        (selected_mzml_paths, selected_idxml_paths, base_directory)
    """
    if header_text:
        st.header(header_text)
    else:
        st.header("📁 Step 2: Select Files")
    
    server_type = server_info.get('type', 'current')
    
    # LOCAL UPLOAD MODE
    if server_type == 'local' and allow_upload:
        st.caption("Upload files from your local computer")
        
        max_size_gb = st.session_state.settings.get("performance", {}).get("max_file_size_gb", 2) if 'settings' in st.session_state else 2
        
        selected_mzml = []
        selected_idxml = []
        base_dir = None
        
        if show_mzml and show_idxml:
            col1, col2 = st.columns(2)
        elif show_mzml:
            col1 = st.container()
            col2 = None
        else:
            col1 = None
            col2 = st.container()
        
        if show_mzml and col1:
            with col1:
                st.subheader("📊 mzML Files")
                uploaded_mzml = st.file_uploader(
                    "Select mzML files",
                    type=["mzml", "mzML"],
                    accept_multiple_files=True,
                    key=f"{session_key_prefix}_mzml_upload",
                    help=f"MS/MS spectra files (up to {max_size_gb} GB each)"
                )
                
                if uploaded_mzml:
                    # Store uploaded files (would need to save to temp in actual processing)
                    selected_mzml = uploaded_mzml  # Returns file objects
                    total_size_mb = sum(f.size for f in uploaded_mzml) / 1024 / 1024
                    st.success(f"✅ {len(uploaded_mzml)} file(s) ({total_size_mb:.1f} MB)")
        
        if show_idxml and col2:
            with col2:
                st.subheader("🔍 idXML Files")
                uploaded_idxml = st.file_uploader(
                    "Select identification files",
                    type=["idxml", "idXML", "mzid", "mzidentml"],
                    accept_multiple_files=True,
                    key=f"{session_key_prefix}_idxml_upload",
                    help=f"Peptide IDs (up to {max_size_gb} GB each)"
                )
                
                if uploaded_idxml:
                    selected_idxml = uploaded_idxml  # Returns file objects
                    total_size_mb = sum(f.size for f in uploaded_idxml) / 1024 / 1024
                    st.success(f"✅ {len(uploaded_idxml)} file(s) ({total_size_mb:.1f} MB)")
        
        return selected_mzml, selected_idxml, base_dir
    
    # SERVER BROWSER MODE (current or ssh servers)
    else:
        st.caption(f"Browse and select files on {server_info.get('display', 'server')}")
        
        # Use the unified file selector (which includes directory browser)
        selected_mzml, selected_idxml, current_dir = render_mzml_idxml_file_selector(
            session_key_prefix=session_key_prefix,
            require_both=require_both,
            show_mzml=show_mzml,
            show_idxml=show_idxml,
            header_text=None,  # We already have a header
            description=None
        )
        
        return selected_mzml, selected_idxml, current_dir


def render_complete_data_input_section(
    session_key_prefix: str = "data_input",
    require_both: bool = False,
    show_mzml: bool = True,
    show_idxml: bool = True,
    allow_upload: bool = True,
    main_header: Optional[str] = None
) -> Tuple[List[str], List[str], str, Dict[str, Any]]:
    """
    Complete data input workflow: Server selection + File selection.
    
    This combines both steps into a unified interface that looks identical
    across all workspaces. Use this for a complete, consistent data input flow.
    
    Args:
        session_key_prefix: Prefix for all session state keys
        require_both: Whether both mzML and idXML are required
        show_mzml: Show mzML file selection
        show_idxml: Show idXML file selection
        allow_upload: Allow file uploads for local mode
        main_header: Main section header
    
    Returns:
        (selected_mzml, selected_idxml, base_directory, server_info)
    """
    if main_header:
        st.header(main_header)
    else:
        st.header("📁 Data Input")
    
    st.markdown("---")
    
    # STEP 1: Server Selection (independent)
    server_info = render_server_selection_section(
        session_key=f"{session_key_prefix}_server",
        header_text="🌐 Step 1: Select Data Source"
    )
    
    st.markdown("---")
    
    # STEP 2: File Selection (independent, based on server)
    selected_mzml, selected_idxml, base_dir = render_file_selection_section(
        server_info=server_info,
        session_key_prefix=f"{session_key_prefix}_files",
        require_both=require_both,
        show_mzml=show_mzml,
        show_idxml=show_idxml,
        header_text="📁 Step 2: Select Files",
        allow_upload=allow_upload
    )
    
    return selected_mzml, selected_idxml, base_dir, server_info
