"""Shared data input UI components for mzML and idXML file selection.

This module provides consistent file selection interfaces used by:
- spectrum_analysis.py: Select files for spectral comparison
- dataset_filter.py: Select files for filtering operations

Usage:
    from data_input_components import render_mzml_idxml_file_selector
    
    mzml_files, idxml_files, base_dir = render_mzml_idxml_file_selector()
"""

import streamlit as st
from pathlib import Path
from typing import Tuple, List, Optional

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
