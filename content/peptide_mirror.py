"""
Peptide Mirror Analysis

Interactive page for comparing query spectra against reference spectra.
"""
import streamlit as st
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple
import json
import logging

from src.common.common import page_setup

params = page_setup()

logger = logging.getLogger(__name__)

# State key prefix for this workspace
STATE_PREFIX = "peptide_mirror"


def format_time(seconds: float) -> str:
    """Format seconds as human-readable time."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"

# ============================================================================
# SETTINGS - loaded by page_setup() into st.session_state.settings
# ============================================================================

def get_settings_value(settings: dict, *keys, default=None):
    """Safely get nested settings value."""
    value = settings
    for key in keys:
        if isinstance(value, dict):
            value = value.get(key)
        else:
            return default
    return value if value is not None else default


def _get_settings() -> dict:
    """Get settings from session state (loaded by page_setup)."""
    return st.session_state.get("settings", {})


def get_state_key(key: str) -> str:
    """Get namespaced session state key."""
    return f"{STATE_PREFIX}.{key}"


def init_session_state():
    """Initialize session state with defaults from settings.json."""
    # Get defaults from settings.json (same values used by main.py)
    ion_types_default = get_settings_value(_get_settings(), "ion_types", "default", default=["b", "y"])
    mz_tolerance_default = get_settings_value(_get_settings(), "mz_tolerance", "default", default=0.02)
    export_dir_default = get_settings_value(_get_settings(), "export", "default_directory", default="./exports")
    export_formats_default = get_settings_value(_get_settings(), "export", "formats", default=["png"])
    project_name_default = get_settings_value(_get_settings(), "export", "project_name", default="MSpecViz_analysis")
    data_dir_default = get_settings_value(_get_settings(), "data_directories", "default", default="/mnt/volume")

    # Annotation settings
    add_losses_default = get_settings_value(_get_settings(), "annotation", "add_losses", default=False)
    add_first_prefix_default = get_settings_value(_get_settings(), "annotation", "add_first_prefix_ion", default=True)
    remove_precursor_default = get_settings_value(_get_settings(), "annotation", "remove_precursor_peaks", default=False)

    workers_default = get_settings_value(_get_settings(), "performance", "workers", default=4)

    defaults = {
        'idxml_dir': data_dir_default,
        'mzml_dir': data_dir_default,
        'selected_idxmls': [],
        'selected_mzmls': [],
        'target_peptides': '',
        'reference_substring': '',
        'ion_types': ion_types_default,
        'mz_tolerance': mz_tolerance_default,
        'add_losses': add_losses_default,
        'add_first_prefix_ion': add_first_prefix_default,
        'remove_precursor_peaks': remove_precursor_default,
        'export_dir': export_dir_default,
        'export_formats': export_formats_default,
        'project_name': project_name_default,
        'workers': workers_default,
        'parquet_file': None,
        'search_results': None,
        'streaming_results': None,
        'current_step': 1,
    }

    for key, default in defaults.items():
        state_key = get_state_key(key)
        if state_key not in st.session_state:
            st.session_state[state_key] = default


def render_file_selection():
    """Render file selection step."""
    st.header("Step 1: File Selection")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(
            key_prefix=f"{STATE_PREFIX}_idxml",
            extensions=['.idXML', '.idxml'],
            state_key='selected_idxmls'
        )

    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(
            key_prefix=f"{STATE_PREFIX}_mzml",
            extensions=['.mzML', '.mzml'],
            state_key='selected_mzmls'
        )

    # Parquet shortcut
    st.markdown("---")
    st.subheader("Pre-converted Parquet (Optional)")
    st.caption("If you have a previously converted Parquet file, select it here to skip search and load spectra instantly.")

    parquet_col1, parquet_col2 = st.columns([3, 1])
    with parquet_col1:
        parquet_path = st.text_input(
            "Parquet file path",
            value=st.session_state.get(get_state_key('parquet_file'), '') or '',
            key=f"{STATE_PREFIX}_parquet_input",
            placeholder="/path/to/converted_spectra.parquet"
        )
    with parquet_col2:
        st.write("")  # spacer
        if parquet_path and Path(parquet_path).exists() and parquet_path.endswith('.parquet'):
            st.session_state[get_state_key('parquet_file')] = parquet_path
            st.success("Parquet file found")
        elif parquet_path:
            st.session_state[get_state_key('parquet_file')] = None
            if not Path(parquet_path).exists():
                st.warning("File not found")
            elif not parquet_path.endswith('.parquet'):
                st.warning("Not a .parquet file")
        else:
            st.session_state[get_state_key('parquet_file')] = None

    # Show selection summary
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    has_parquet = bool(st.session_state.get(get_state_key('parquet_file')))

    if has_parquet:
        st.success("Parquet file selected — skip to **Process** tab (Steps 2-5 are optional)")
        return True
    elif idxmls and mzmls:
        st.success(f"Selected {len(idxmls)} idXML file(s) and {len(mzmls)} mzML file(s)")
        return True
    else:
        if not idxmls:
            st.warning("Please select at least one idXML file")
        if not mzmls:
            st.warning("Please select at least one mzML file")
        return False


def render_directory_browser(key_prefix: str, extensions: list, state_key: str):
    """Render a directory browser with file selection."""
    # State keys
    dir_key = f"{key_prefix}_dir"
    nav_key = f"{key_prefix}_nav_to"  # Pending navigation target
    select_all_key = f"{key_prefix}_select_all"

    # Initialize directory state
    if dir_key not in st.session_state:
        st.session_state[dir_key] = "/mnt/volume"

    # Handle pending navigation (set by button callbacks)
    if nav_key in st.session_state and st.session_state[nav_key]:
        st.session_state[dir_key] = st.session_state[nav_key]
        st.session_state[nav_key] = None
        # Reset select all when changing directory
        if select_all_key in st.session_state:
            del st.session_state[select_all_key]

    current_dir = st.session_state[dir_key]

    # Path input with up button
    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path",
            value=current_dir,
            key=f"{key_prefix}_path_input_{hash(current_dir) % 10000}",  # Dynamic key forces refresh
            label_visibility="collapsed",
            placeholder="Enter path..."
        )
        # Handle manual path entry
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[dir_key] = new_path
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("⬆️", key=f"{key_prefix}_up_btn", help="Parent directory", disabled=(parent == current_dir)):
            st.session_state[nav_key] = parent
            st.rerun()

    # Show clickable subdirectories
    try:
        subdirs = sorted([d for d in Path(current_dir).iterdir() if d.is_dir()])
        if subdirs:
            st.caption("Click folder to navigate:")
            # Display folders in a grid
            cols_per_row = 3
            max_folders = 15
            for i in range(0, min(len(subdirs), max_folders), cols_per_row):
                cols = st.columns(cols_per_row)
                for j, col in enumerate(cols):
                    idx = i + j
                    if idx < len(subdirs) and idx < max_folders:
                        folder = subdirs[idx]
                        folder_path = str(folder)
                        with col:
                            if st.button(
                                f"📁 {folder.name}",
                                key=f"{key_prefix}_folder_{idx}",
                                width="stretch"
                            ):
                                st.session_state[nav_key] = folder_path
                                st.rerun()
            if len(subdirs) > max_folders:
                st.caption(f"... and {len(subdirs) - max_folders} more folders")
    except PermissionError:
        st.warning("Permission denied")
    except Exception as e:
        st.error(f"Cannot access: {e}")

    st.markdown("---")

    # Get matching files
    try:
        all_files = []
        for ext in extensions:
            all_files.extend(Path(current_dir).glob(f"*{ext}"))
        all_files = sorted(all_files, key=lambda x: x.name)

        # Filter input
        filter_text = st.text_input(
            "Filter files",
            key=f"{key_prefix}_filter",
            placeholder="Type to filter filenames..."
        )

        # Apply filter
        if filter_text:
            filter_lower = filter_text.lower()
            filtered_files = [f for f in all_files if filter_lower in f.name.lower()]
        else:
            filtered_files = all_files

        if filtered_files:
            file_names = [f.name for f in filtered_files]

            # Get ALL currently selected files (from all directories)
            all_selected_paths = st.session_state.get(get_state_key(state_key), [])

            # Separate selections: current directory vs other directories
            current_dir_selections = [p for p in all_selected_paths if Path(p).parent == Path(current_dir)]
            other_dir_selections = [p for p in all_selected_paths if Path(p).parent != Path(current_dir)]

            # Get names of files selected in current directory (for multiselect default)
            current_dir_names = [Path(p).name for p in current_dir_selections if Path(p).exists()]
            valid_defaults = [n for n in current_dir_names if n in file_names]

            # Check if all files in current directory are selected
            all_in_dir_selected = set(valid_defaults) == set(file_names) and len(valid_defaults) > 0

            # Select All checkbox (for current directory only)
            select_all = st.checkbox(
                f"Select all {len(filtered_files)} files" + (f" (filtered from {len(all_files)})" if filter_text else ""),
                value=all_in_dir_selected,
                key=select_all_key
            )

            if select_all:
                # All files in current directory selected
                selected_names = file_names
            else:
                # Manual selection via multiselect
                selected_names = st.multiselect(
                    "Select files from this folder:",
                    options=file_names,
                    default=valid_defaults if not all_in_dir_selected else [],
                    key=f"{key_prefix}_files"
                )

            # Build new selection: other directories + current directory selections
            # Deduplicate exact same paths (preserves order)
            current_dir_selected_paths = [str(f) for f in filtered_files if f.name in selected_names]
            merged_selections = list(dict.fromkeys(other_dir_selections + current_dir_selected_paths))
            st.session_state[get_state_key(state_key)] = merged_selections

            if current_dir_selected_paths:
                st.success(f"✅ {len(current_dir_selected_paths)} file(s) selected from this folder")
        else:
            if filter_text:
                st.info(f"No files match filter '{filter_text}'")
            else:
                st.info(f"No {extensions[0]} files in this directory")

    except Exception as e:
        st.error(f"Error listing files: {e}")

    # Display accumulated selections from all folders
    all_selected = st.session_state.get(get_state_key(state_key), [])
    if all_selected:
        # Group by directory
        grouped = defaultdict(list)
        for path in all_selected:
            p = Path(path)
            if p.exists():
                grouped[str(p.parent)].append(p.name)

        if grouped:
            st.markdown("---")
            with st.expander(f"📂 All selected files ({len(all_selected)} total)", expanded=False):
                # Clear all button
                if st.button("🗑️ Clear all selections", key=f"{key_prefix}_clear_all"):
                    st.session_state[get_state_key(state_key)] = []
                    st.rerun()

                # Show grouped files with remove buttons
                remove_idx = 0
                for folder, files in sorted(grouped.items()):
                    st.caption(f"**{folder}/** ({len(files)} files)")
                    for fname in sorted(files):
                        col1, col2 = st.columns([6, 1])
                        with col1:
                            st.text(f"  {fname}")
                        with col2:
                            full_path = str(Path(folder) / fname)
                            if st.button("❌", key=f"{key_prefix}_remove_{remove_idx}", help=f"Remove {fname}"):
                                updated = [p for p in all_selected if p != full_path]
                                st.session_state[get_state_key(state_key)] = updated
                                st.rerun()
                        remove_idx += 1


def render_peptide_input():
    """Render peptide input step."""
    st.header("Step 2: Target Peptides")

    col1, col2 = st.columns([2, 1])

    with col1:
        peptide_text = st.text_area(
            "Enter peptide sequences (one per line)",
            value=st.session_state.get(get_state_key('target_peptides'), ''),
            height=150,
            key=f"{STATE_PREFIX}_peptide_input",
            placeholder="PEPTIDESEQ\nANOTHERPEPTIDE\n..."
        )
        st.session_state[get_state_key('target_peptides')] = peptide_text

    with col2:
        st.markdown("**Tips:**")
        st.markdown("""
        - Enter base sequences (no modifications)
        - One peptide per line
        - Modified variants will be found automatically
        - Case insensitive
        """)

        # Count entered peptides
        peptides = [p.strip() for p in peptide_text.split('\n') if p.strip()]
        if peptides:
            st.success(f"{len(peptides)} peptide(s) entered")

            # Preview button
            if st.button("🔍 Preview in idXML", key=f"{STATE_PREFIX}_preview_btn"):
                with st.spinner("Scanning idXML files..."):
                    preview_peptides_in_idxmls()

    return len(peptides) > 0


def preview_peptides_in_idxmls():
    """Preview which peptides are found in selected idXML files."""
    from src.peptide_search import get_unique_peptides_from_idxmls

    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    if not idxmls:
        st.warning("No idXML files selected")
        return

    # Get sample of peptides from files
    workers = st.session_state.get(get_state_key('workers'), 4)
    sample_peptides = get_unique_peptides_from_idxmls(idxmls, limit=50, workers=workers)

    if sample_peptides:
        st.info(f"Sample peptides from idXML files (first 50):")
        # Display in columns
        cols = st.columns(3)
        for i, pep in enumerate(sample_peptides[:30]):
            cols[i % 3].caption(pep)


def render_reference_selection():
    """Render reference file selection step."""
    st.header("Step 3: Reference Identification")

    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    if not mzmls:
        st.warning("Please select mzML files first")
        return False

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")

    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=st.session_state.get(get_state_key('reference_substring'), ''),
        key=f"{STATE_PREFIX}_ref_input",
        placeholder="e.g., Synpep"
    )
    st.session_state[get_state_key('reference_substring')] = ref_substring

    # Show categorization preview
    if ref_substring:
        reference_files = []
        query_files = []

        for mzml in mzmls:
            fname = Path(mzml).name
            if ref_substring.lower() in fname.lower():
                reference_files.append(fname)
            else:
                query_files.append(fname)

        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Reference files ({len(reference_files)})**")
            for f in reference_files[:5]:
                st.caption(f"📊 {f}")
            if len(reference_files) > 5:
                st.caption(f"... and {len(reference_files) - 5} more")

        with col2:
            st.markdown(f"**Query files ({len(query_files)})**")
            for f in query_files[:5]:
                st.caption(f"📄 {f}")
            if len(query_files) > 5:
                st.caption(f"... and {len(query_files) - 5} more")

        if reference_files and query_files:
            st.success("Reference/Query separation configured")
            return True
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")

    else:
        st.info("No reference identifier — each spectrum will be plotted individually with annotations.")
        return True

    return False


def render_settings():
    """Render analysis settings."""
    st.header("Step 4: Settings")

    # Get available ion types from settings
    available_ions = get_settings_value(_get_settings(), "ion_types", "available",
                                        default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Ion Types")
        current_ions = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])

        selected_ions = st.multiselect(
            "Select ion types to annotate",
            options=available_ions,
            default=[i for i in current_ions if i in available_ions],
            key=f"{STATE_PREFIX}_ion_select"
        )
        st.session_state[get_state_key('ion_types')] = selected_ions

    with col2:
        st.subheader("m/z Tolerance")
        current_tol = st.session_state.get(get_state_key('mz_tolerance'), 0.02)

        tolerance = st.number_input(
            "Tolerance (Da)",
            min_value=0.001,
            max_value=1.0,
            value=float(current_tol),
            step=0.005,
            format="%.3f",
            key=f"{STATE_PREFIX}_tol_input"
        )
        st.session_state[get_state_key('mz_tolerance')] = tolerance

    # Additional annotation settings
    st.markdown("---")
    st.subheader("Annotation Settings")

    ann_col1, ann_col2, ann_col3 = st.columns(3)

    with ann_col1:
        add_losses = st.checkbox(
            "Add neutral losses (H₂O, NH₃)",
            value=st.session_state.get(get_state_key('add_losses'), False),
            key=f"{STATE_PREFIX}_add_losses",
            help="Include neutral loss peaks in theoretical spectrum"
        )
        st.session_state[get_state_key('add_losses')] = add_losses

    with ann_col2:
        add_first_prefix = st.checkbox(
            "Add first prefix ion (a₁, b₁)",
            value=st.session_state.get(get_state_key('add_first_prefix_ion'), True),
            key=f"{STATE_PREFIX}_add_first_prefix",
            help="Include first prefix ions in theoretical spectrum"
        )
        st.session_state[get_state_key('add_first_prefix_ion')] = add_first_prefix

    with ann_col3:
        remove_precursor = st.checkbox(
            "Remove precursor peaks from plot",
            value=st.session_state.get(get_state_key('remove_precursor_peaks'), False),
            key=f"{STATE_PREFIX}_remove_precursor",
            help="Remove all precursor ions (at any charge state) from the final plot. "
                 "Intensities will be renormalized after removal."
        )
        st.session_state[get_state_key('remove_precursor_peaks')] = remove_precursor

    # Intensity mode
    st.markdown("---")
    st.subheader("Intensity Display")
    intensity_mode = st.radio(
        "Intensity mode",
        options=["Relative", "Raw"],
        index=0 if st.session_state.get(get_state_key('intensity_mode'), 'Relative') == 'Relative' else 1,
        horizontal=True,
        key=f"{STATE_PREFIX}_intensity_mode",
        help="Relative: each spectrum normalized to its base peak (0–1). "
             "Raw: original intensity values from the mzML file."
    )
    st.session_state[get_state_key('intensity_mode')] = intensity_mode

    # Info about precursor behavior
    if remove_precursor:
        st.caption("ℹ️ Precursor peaks will be removed from plots. Relative intensities will be recalculated.")
    elif "precursor" in st.session_state.get(get_state_key('ion_types'), []):
        st.caption("ℹ️ Precursor peaks will be shown and colored.")
    else:
        st.caption("ℹ️ Precursor peaks will be shown in grey (unannotated).")

    # Performance settings
    st.markdown("---")
    st.subheader("Performance")
    workers = st.slider(
        "Parallel workers",
        min_value=1,
        max_value=32,
        value=st.session_state.get(get_state_key('workers'), 4),
        key=f"{STATE_PREFIX}_workers",
        help="Number of parallel processes for plot generation. "
             "Higher values speed up export but use more memory (~30 MB per worker)."
    )
    st.session_state[get_state_key('workers')] = workers

    return True


def render_search():
    """Render peptide search step."""
    st.header("Step 5: Search Peptides")

    # Check prerequisites - read from both possible state keys for robustness
    # Files may be stored directly or via the widget
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])

    # Peptides: try the synced key first, then the widget key directly
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    if not peptide_text:
        # Fallback: read directly from text_area widget key
        peptide_text = st.session_state.get(f"{STATE_PREFIX}_peptide_input", '')
    peptides = [p.strip().upper() for p in peptide_text.split('\n') if p.strip()]

    if not idxmls:
        st.warning("Please select idXML files first (Files tab)")
        return

    if not peptides:
        st.warning("Please enter target peptides first (Peptides tab)")
        return

    # Show what will be searched
    st.info(f"Ready to search {len(peptides)} peptide(s) in {len(idxmls)} idXML file(s)")

    if st.button("🔍 Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(idxmls, peptides)

    # Display search results
    search_results = st.session_state.get(get_state_key('search_results'))
    if search_results:
        st.markdown("---")
        st.subheader("Search Results")
        display_search_results(search_results)

        total_matches = sum(len(vg.matches) for vg in search_results.values())
        st.success(f"Found {total_matches} PSMs across {len(search_results)} peptide variants. Go to **Analysis** tab to load spectra.")

        # Convert to Parquet button
        mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
        ref_sub = st.session_state.get(get_state_key('reference_substring'), '')
        if mzmls:
            st.markdown("---")
            st.subheader("Convert to Parquet (Optional)")
            st.caption("Pre-convert matched spectra to Parquet for instant loading on future runs.")

            parquet_out_col1, parquet_out_col2 = st.columns([3, 1])
            with parquet_out_col1:
                default_parquet_dir = st.session_state.get(get_state_key('export_dir'), './exports')
                parquet_output = st.text_input(
                    "Parquet output path",
                    value=f"{default_parquet_dir}/spectra_cache.parquet",
                    key=f"{STATE_PREFIX}_parquet_output",
                )
            with parquet_out_col2:
                st.write("")
                if st.button("Convert to Parquet", key=f"{STATE_PREFIX}_convert_parquet_btn"):
                    _run_parquet_conversion(search_results, mzmls, ref_sub, parquet_output)
    else:
        st.caption("Search results will appear here after running the search.")


def _run_parquet_conversion(search_results, mzmls, ref_substring, output_path):
    """Run inline Parquet conversion from search results."""
    from tools.convert_mzml_to_parquet import convert_mzml_to_parquet

    # Collect unique idXML resolved paths and peptides from search results
    idxml_resolved_paths = set()
    peptides = set()
    for vg in search_results.values():
        peptides.add(vg.base_sequence)
        for m in vg.matches:
            if m.idxml_resolved_path:
                idxml_resolved_paths.add(m.idxml_resolved_path)

    # Determine mzML directory (use common parent of selected mzML files)
    mzml_dirs = {str(Path(m).parent) for m in mzmls}
    if len(mzml_dirs) == 1:
        mzml_dir = mzml_dirs.pop()
    else:
        # Use first directory as primary
        mzml_dir = str(Path(mzmls[0]).parent)

    # Use resolved paths directly (no name-based lookup needed)
    idxml_paths = list(idxml_resolved_paths)

    # Fallback to name-based lookup if resolved paths not available
    if not idxml_paths:
        idxml_names = set()
        for vg in search_results.values():
            for m in vg.matches:
                if m.idxml_file:
                    idxml_names.add(m.idxml_file)
        selected_idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
        for idxml_name in idxml_names:
            for full_path in selected_idxmls:
                if Path(full_path).name == idxml_name:
                    idxml_paths.append(full_path)
                    break

    if not idxml_paths:
        st.error("Could not find full paths for idXML files")
        return

    progress_bar = st.progress(0, "Converting...")
    status = st.empty()

    def on_progress(msg):
        status.caption(msg)

    def on_error(msg):
        st.warning(msg)

    result = convert_mzml_to_parquet(
        idxml_paths=idxml_paths,
        mzml_dir=mzml_dir,
        target_peptides=list(peptides),
        output_path=output_path,
        reference_substring=ref_substring,
        progress_callback=on_progress,
        error_callback=on_error,
    )

    progress_bar.progress(1.0)

    if result:
        st.session_state[get_state_key('parquet_file')] = result
        st.success(f"Converted to: {result}")
    else:
        st.error("Conversion failed")


def _show_cache_status(export_dir: str, parquet_file: str):
    """Show whether a Parquet cache exists for current inputs."""
    if parquet_file:
        # User manually specified a Parquet file — cache is irrelevant
        return

    settings = _get_settings()
    cache_enabled = get_settings_value(settings, "cache", "enabled", default=True)
    if not cache_enabled:
        return

    cache_subdir = get_settings_value(settings, "cache", "subdirectory", default=".mspecviz_cache")

    # Gather inputs for fingerprint
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    peptides = [p.strip().upper() for p in peptide_text.split('\n') if p.strip()]
    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')

    if not idxmls or not mzmls or not peptides:
        return

    from src.parquet_cache import compute_cache_fingerprint, find_existing_cache, is_cache_complete
    fingerprint = compute_cache_fingerprint(idxmls, mzmls, peptides, ref_substring)
    cache_dir = str(Path(export_dir) / cache_subdir)
    cached = find_existing_cache(cache_dir, fingerprint)

    if cached:
        st.success(f"Parquet cache found — loading will be instant ({Path(cached).name})")
    else:
        # Check for partial cache
        partial = find_existing_cache(cache_dir, fingerprint, allow_partial=True)
        if partial:
            from src.parquet_cache import get_parquet_row_counts
            refs, queries = get_parquet_row_counts(partial)
            st.info(f"Partial cache found ({refs} refs, {queries} queries). "
                    f"Resume will load remaining spectra from mzML.")
        else:
            st.caption("No cache yet. First run creates a Parquet cache for instant loading next time.")


def render_analysis():
    """Render spectrum loading and analysis step."""
    st.header("Step 6: Process & Export")

    # Check prerequisites — Parquet path OR (search_results + mzmls)
    parquet_file = st.session_state.get(get_state_key('parquet_file'))
    search_results = st.session_state.get(get_state_key('search_results'))
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')
    if not ref_substring:
        ref_substring = st.session_state.get(f"{STATE_PREFIX}_ref_input", '')

    if not parquet_file and not search_results:
        st.warning("Please run peptide search first (Search tab) or select a Parquet file (Files tab)")
        return

    if not parquet_file and not mzmls:
        st.warning("Please select mzML files first (Files tab)")
        return

    # Show mode indicator
    if parquet_file:
        st.info(f"Loading from Parquet: {Path(parquet_file).name}")
    else:
        total_psms = sum(len(vg.matches) for vg in search_results.values())

    streaming_results = st.session_state.get(get_state_key('streaming_results'))

    # Show status as caption (less prominent)
    if streaming_results:
        st.caption(f"✅ Processed {streaming_results.get('query_spectra', 0)} spectra, "
                   f"exported {streaming_results.get('plots_generated', 0)} plots")
    elif not parquet_file:
        st.caption(f"Ready to process {total_psms} peptide-spectrum matches")

    st.markdown("---")

    # Export settings
    with st.expander("📁 Export Settings", expanded=True):
        export_col1, export_col2 = st.columns(2)

        with export_col1:
            project_name = st.text_input(
                "Project Name",
                value=st.session_state.get(get_state_key('project_name'), 'MSpecViz_analysis'),
                key=f"{STATE_PREFIX}_project_name_input"
            )
            st.session_state[get_state_key('project_name')] = project_name

            export_dir = st.text_input(
                "Export Directory",
                value=st.session_state.get(get_state_key('export_dir'), './exports'),
                key=f"{STATE_PREFIX}_export_dir_input"
            )
            st.session_state[get_state_key('export_dir')] = export_dir

        with export_col2:
            export_formats = st.multiselect(
                "Export Formats",
                options=['png'],
                default=st.session_state.get(get_state_key('export_formats'), ['png']),
                key=f"{STATE_PREFIX}_export_formats_input"
            )
            st.session_state[get_state_key('export_formats')] = export_formats

            plot_unpaired = st.checkbox(
                "Plot unpaired spectra",
                value=st.session_state.get(get_state_key('plot_unpaired'), False),
                key=f"{STATE_PREFIX}_plot_unpaired_input",
                help="When enabled, spectra without a matching reference (or query) are plotted individually. "
                     "When disabled, only mirror plots with both query and reference are exported."
            )
            st.session_state[get_state_key('plot_unpaired')] = plot_unpaired

    # Cache status indicator
    _show_cache_status(export_dir, parquet_file)

    # Process button
    if st.button("🚀 Process & Export", key=f"{STATE_PREFIX}_load_btn", type="primary",
                 help="Load spectra, generate plots, and export"):
        run_streaming_processing()

    # Display results
    if streaming_results:
        st.markdown("---")
        st.subheader("Results")
        display_streaming_results(streaming_results)
    else:
        st.caption("Click 'Process & Export' to load spectra, generate plots, and export.")


def run_peptide_search(idxmls: list, peptides: list):
    """Execute peptide search (with checkpoint resume support)."""
    from src.peptide_search import search_peptides_in_idxmls

    # Get selected mzML paths for disambiguation
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    # Check for search checkpoint before running the full search
    settings = _get_settings()
    cache_enabled = get_settings_value(settings, "cache", "enabled", default=True)
    if cache_enabled and mzmls:
        cache_subdir = get_settings_value(settings, "cache", "subdirectory", default=".mspecviz_cache")
        export_dir_ck = st.session_state.get(get_state_key('export_dir'), './exports')
        ref_sub = st.session_state.get(get_state_key('reference_substring'), '')
        peptide_text_ck = st.session_state.get(get_state_key('target_peptides'), '')
        peptides_ck = [p.strip().upper() for p in peptide_text_ck.split('\n') if p.strip()]

        from src.parquet_cache import compute_cache_fingerprint, load_search_checkpoint
        ck_fingerprint = compute_cache_fingerprint(idxmls, mzmls, peptides_ck, ref_sub)
        ck_cache_dir = str(Path(export_dir_ck) / cache_subdir)
        checkpoint_results = load_search_checkpoint(ck_cache_dir, ck_fingerprint)
        if checkpoint_results:
            total_psms = sum(len(vg.matches) for vg in checkpoint_results.values())
            st.session_state[get_state_key('search_results')] = checkpoint_results
            st.success(f"Search checkpoint loaded: {total_psms} PSMs across "
                       f"{len(checkpoint_results)} variants (skipped idXML search)")
            return

    error_messages = []

    import time as _time

    with st.status("Searching peptides...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        _search_start = _time.time()

        def progress_callback(current, total, message):
            if total > 0:
                frac = current / total
                progress_bar.progress(frac, message)
                elapsed = _time.time() - _search_start
                if frac > 0 and current < total:
                    remaining = (elapsed / frac) - elapsed
                    status_text.caption(
                        f"{current}/{total} files  "
                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)"
                    )
                elif current >= total:
                    status_text.caption(f"Done in {format_time(elapsed)}")

        def error_callback(error_msg: str):
            """Handle errors immediately."""
            error_messages.append(error_msg)

        try:
            workers = st.session_state.get(get_state_key('workers'), 4)
            results = search_peptides_in_idxmls(
                idxmls, peptides,
                mzml_paths=mzmls,
                progress_callback=progress_callback,
                error_callback=error_callback,
                workers=workers,
            )

            st.session_state[get_state_key('search_results')] = results
            progress_bar.progress(1.0, "Search complete!")

            if results:
                total_matches = sum(len(vg.matches) for vg in results.values())

                # Save search checkpoint immediately
                settings = _get_settings()
                cache_enabled = get_settings_value(settings, "cache", "enabled", default=True)
                if cache_enabled:
                    cache_subdir = get_settings_value(settings, "cache", "subdirectory", default=".mspecviz_cache")
                    export_dir_ck = st.session_state.get(get_state_key('export_dir'), './exports')
                    ref_sub = st.session_state.get(get_state_key('reference_substring'), '')
                    peptide_text_ck = st.session_state.get(get_state_key('target_peptides'), '')
                    peptides_ck = [p.strip().upper() for p in peptide_text_ck.split('\n') if p.strip()]

                    from src.parquet_cache import compute_cache_fingerprint, save_search_checkpoint
                    ck_fingerprint = compute_cache_fingerprint(idxmls, mzmls, peptides_ck, ref_sub)
                    ck_cache_dir = str(Path(export_dir_ck) / cache_subdir)
                    ck_saved = save_search_checkpoint(ck_cache_dir, ck_fingerprint, results)
                    if ck_saved:
                        st.caption(f"Search checkpoint saved: {Path(ck_saved).name}")

                # Count resolved vs unresolved mzML mappings
                resolved = sum(1 for vg in results.values() for m in vg.matches if m.mzml_filename)
                unresolved = total_matches - resolved

                status.update(label=f"Found {total_matches} PSMs across {len(results)} variants", state="complete")

                if unresolved > 0:
                    st.warning(f"{unresolved} PSMs could not be mapped to mzML files (missing id_merge_index)")

                # Show error summary
                if error_messages:
                    with st.expander(f"{len(error_messages)} warnings during search", expanded=False):
                        for msg in error_messages[:10]:
                            st.caption(msg)
                        if len(error_messages) > 10:
                            st.caption(f"... and {len(error_messages) - 10} more")
            else:
                status.update(label="No matches found for the specified peptides", state="error")

        except Exception as e:
            status.update(label=f"Search failed: {e}", state="error")
            logger.exception("Peptide search error")


def run_streaming_processing():
    """
    Run streaming spectrum loading and processing.

    Supports two paths:
    - Parquet fast path: load from pre-converted .parquet file
    - mzML path: stream from mzML files with optimized index building
    """
    from src.spectrum_loader import load_spectra_streaming, load_spectra_from_parquet, StreamingResult
    from src.peak_annotation import prepare_spectrum_for_plot
    from src.export import sanitize_filename, shorten_stem, create_export_directory
    from src.spectrum_metrics import shorten_filenames
    from src.parallel_processor import ComparisonTask, run_parallel_comparisons

    # Get state
    parquet_file = st.session_state.get(get_state_key('parquet_file'))
    search_results = st.session_state.get(get_state_key('search_results'))
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')
    if not ref_substring:
        ref_substring = st.session_state.get(f"{STATE_PREFIX}_ref_input", '')

    ion_types = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])
    mz_tolerance = st.session_state.get(get_state_key('mz_tolerance'), 0.02)
    add_losses = st.session_state.get(get_state_key('add_losses'), False)
    add_first_prefix_ion = st.session_state.get(get_state_key('add_first_prefix_ion'), True)
    remove_precursor_peaks = st.session_state.get(get_state_key('remove_precursor_peaks'), False)
    intensity_mode = st.session_state.get(get_state_key('intensity_mode'), 'Relative')
    normalize = intensity_mode == 'Relative'
    export_dir = st.session_state.get(get_state_key('export_dir'), './exports')
    export_formats = st.session_state.get(get_state_key('export_formats'), ['png'])
    project_name = st.session_state.get(get_state_key('project_name'), 'MSpecViz_analysis')
    plot_unpaired = st.session_state.get(get_state_key('plot_unpaired'), False)

    # =========================================================================
    # Auto-cache detection: check for existing Parquet cache before loading
    # =========================================================================
    loaded_from_cache = False
    cache_dir = None
    cache_fingerprint = None

    settings = _get_settings()
    cache_enabled = get_settings_value(settings, "cache", "enabled", default=True)
    cache_subdir = get_settings_value(settings, "cache", "subdirectory", default=".mspecviz_cache")

    partial_cache_path = None

    if not parquet_file and cache_enabled:
        # Only auto-detect when user hasn't manually specified a Parquet file
        idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
        peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
        peptides_list = [p.strip().upper() for p in peptide_text.split('\n') if p.strip()]

        if idxmls and mzmls and peptides_list:
            from src.parquet_cache import compute_cache_fingerprint, find_existing_cache
            cache_fingerprint = compute_cache_fingerprint(idxmls, mzmls, peptides_list, ref_substring)
            cache_dir = str(Path(export_dir) / cache_subdir)

            # Check for complete cache first
            cached_path = find_existing_cache(cache_dir, cache_fingerprint)
            if cached_path:
                parquet_file = cached_path
                loaded_from_cache = True
                st.info(f"Cache found! Loading from: {Path(cached_path).name}")
            else:
                # Check for partial (incomplete) cache
                partial_path = find_existing_cache(cache_dir, cache_fingerprint, allow_partial=True)
                if partial_path:
                    from src.parquet_cache import get_parquet_row_counts
                    p_refs, p_queries = get_parquet_row_counts(partial_path)
                    partial_cache_path = partial_path
                    st.info(f"Partial cache found ({p_refs} refs, {p_queries} queries). "
                            f"Resuming — loading remaining spectra from mzML.")

    use_parquet = bool(parquet_file)

    if not use_parquet:
        if not search_results:
            st.warning("No search results to process")
            return

        # Make a copy of search results
        import copy
        variant_groups = copy.deepcopy(search_results)
        mzml_paths = [Path(m) for m in mzmls]
        total_psms = sum(len(vg.matches) for vg in variant_groups.values())
    else:
        variant_groups = None
        mzml_paths = None
        total_psms = 0  # Will be determined from Parquet

    # Build run settings dict (everything that affects output)
    from src.export import find_or_create_export_directory

    dpi = get_settings_value(_get_settings(), "export", "dpi", default=300)
    figsize = list(get_settings_value(_get_settings(), "export", "figsize", default=[12, 8]))

    idxmls_for_settings = st.session_state.get(get_state_key('selected_idxmls'), [])
    peptide_text_s = st.session_state.get(get_state_key('target_peptides'), '')
    peptides_for_settings = sorted(p.strip().upper() for p in peptide_text_s.split('\n') if p.strip())

    run_settings = {
        "idxml_files": sorted(idxmls_for_settings),
        "mzml_files": sorted(mzmls),
        "peptides": peptides_for_settings,
        "reference_substring": ref_substring,
        "ion_types": sorted(ion_types),
        "mz_tolerance": mz_tolerance,
        "add_losses": add_losses,
        "add_first_prefix_ion": add_first_prefix_ion,
        "remove_precursor_peaks": remove_precursor_peaks,
        "intensity_mode": intensity_mode,
        "export_formats": sorted(export_formats),
        "dpi": dpi,
        "figsize": figsize,
        "plot_unpaired": plot_unpaired,
    }

    # Numbered folder: MSpecViz_analysis_001/, _002/, etc. Reuses matching folder.
    root_dir, png_dir, svg_dir = find_or_create_export_directory(
        export_dir, project_name, run_settings
    )

    # Count expected reference vs query PSMs
    if use_parquet:
        # Read row counts from Parquet for accurate progress bars
        from src.parquet_cache import get_parquet_row_counts
        expected_refs, expected_queries = get_parquet_row_counts(parquet_file)
    else:
        # Count refs using resolved paths when available, bare name as fallback
        ref_paths = {str(Path(m).resolve()).lower() for m in mzmls
                     if ref_substring and ref_substring.lower() in Path(m).name.lower()}
        ref_names_lower = {Path(m).name.lower() for m in mzmls
                           if ref_substring and ref_substring.lower() in Path(m).name.lower()}
        expected_refs = sum(
            1 for vg in variant_groups.values() for m in vg.matches
            if (m.mzml_resolved_path and str(Path(m.mzml_resolved_path).resolve()).lower() in ref_paths)
            or (not m.mzml_resolved_path and m.mzml_filename and m.mzml_filename.lower() in ref_names_lower)
        )
        expected_queries = total_psms - expected_refs

    # Progress display - two separate progress bars
    import time
    start_time = time.time()
    ref_start_time = start_time
    query_start_time = None

    # Reference progress section
    refs_done = False
    if expected_refs > 0:
        st.markdown("**Reference Spectra**")
        ref_progress_bar = st.progress(0)
        ref_status = st.empty()
        ref_status.caption(f"Loading {expected_refs} reference spectra...")
    else:
        ref_progress_bar = None
        ref_status = None
        refs_done = True
        query_start_time = start_time
        st.info("No reference substring specified — generating individual annotated spectrum plots.")

    # Query progress section
    st.markdown("**Query Spectra**" if expected_refs > 0 else "**Spectra**")
    query_progress_bar = st.progress(0)
    query_status = st.empty()
    if expected_refs > 0:
        query_status.caption(f"Waiting for references to complete...")
    else:
        query_status.caption(f"Processing {expected_queries} spectra...")

    # Summary metrics (compact)
    metrics_col1, metrics_col2 = st.columns(2)
    plots_metric = metrics_col1.empty()
    errors_metric = metrics_col2.empty()

    def add_log(message: str):
        """Log to terminal (not app UI)."""
        logger.info(message)
        print(message)

    # Results tracking (same structure as main.py)
    results_summary = {
        "reference_spectra": 0,
        "query_spectra": 0,
        "plots_generated": 0,
        "files_exported": {"png": [], "svg": []},
        "errors": [],
        "export_root": str(root_dir),
    }

    # Annotated reference cache (annotate once, reuse)
    annotated_references = {}
    # refs_done initialized in the reference progress section above

    # Accumulate per-comparison rows for summary spreadsheet
    summary_rows: List[dict] = []

    # Accumulate spectra for auto-cache save (only when loading from mzML)
    _cached_queries: List[Tuple[str, object]] = []
    _cached_references: Dict[str, list] = {}

    # Incremental cache save state
    _INCREMENTAL_SAVE_INTERVAL = 500
    _last_incremental_save = [0]
    _cache_complete = [False]

    def _try_save_cache(stage: str, partial: bool = False):
        """Save Parquet cache. Supports incremental partial saves."""
        if _cache_complete[0]:
            return
        if loaded_from_cache or st.session_state.get(get_state_key('parquet_file')):
            return
        if not cache_enabled or not cache_fingerprint:
            return
        has_refs = _cached_references and any(_cached_references.values())
        if not has_refs and not _cached_queries:
            return
        try:
            from src.parquet_cache import save_spectra_to_parquet
            c_path = str(Path(cache_dir) / f"{cache_fingerprint}.parquet")

            idxmls_for_meta = st.session_state.get(get_state_key('selected_idxmls'), [])
            peptide_text_for_meta = st.session_state.get(get_state_key('target_peptides'), '')
            peptides_for_meta = [p.strip().upper() for p in peptide_text_for_meta.split('\n') if p.strip()]

            cache_metadata = {
                "fingerprint": cache_fingerprint,
                "export_root": str(root_dir),
                "idxml_files": idxmls_for_meta,
                "mzml_files": mzmls,
                "target_peptides": peptides_for_meta,
                "reference_substring": ref_substring,
            }

            saved = save_spectra_to_parquet(
                reference_spectra=_cached_references,
                query_spectra=_cached_queries,
                output_path=c_path,
                cache_metadata=cache_metadata,
                partial=partial,
            )
            if saved:
                if not partial:
                    _cache_complete[0] = True
                tag = "partial" if partial else "complete"
                add_log(f"Cache saved [{tag}] ({stage}): {Path(saved).name}")
                if not partial:
                    st.caption(f"Cache saved for instant loading next time: {Path(saved).name}")
        except Exception as cache_err:
            add_log(f"Warning: cache save failed ({stage}): {cache_err}")
            logger.warning(f"Cache save failed ({stage}): {cache_err}")

    # Task collection for parallel processing
    comparison_tasks: List[ComparisonTask] = []
    task_counter = [0]

    def annotate_spectrum(spectrum, is_reference=False):
        """Annotate a spectrum for plotting (same as main.py)."""
        spec_type = "reference" if is_reference else "query"

        # Validate input data
        if spectrum.mz is None or len(spectrum.mz) == 0:
            error_msg = f"Annotation error ({spec_type} {spectrum.peptidoform}): Empty m/z array"
            results_summary["errors"].append(error_msg)
            add_log(f"    ❌ {error_msg}")
            return None

        if spectrum.intensity is None or len(spectrum.intensity) == 0:
            error_msg = f"Annotation error ({spec_type} {spectrum.peptidoform}): Empty intensity array"
            results_summary["errors"].append(error_msg)
            add_log(f"    ❌ {error_msg}")
            return None

        if not spectrum.peptidoform:
            error_msg = f"Annotation error ({spec_type}): Missing peptidoform"
            results_summary["errors"].append(error_msg)
            add_log(f"    ❌ {error_msg}")
            return None

        if spectrum.precursor_charge <= 0:
            error_msg = f"Annotation error ({spec_type} {spectrum.peptidoform}): Invalid charge {spectrum.precursor_charge}"
            results_summary["errors"].append(error_msg)
            add_log(f"    ❌ {error_msg}")
            return None

        try:
            df = prepare_spectrum_for_plot(
                mz=spectrum.mz,
                intensity=spectrum.intensity,
                peptidoform=spectrum.peptidoform,
                charge=spectrum.precursor_charge,
                ion_types=ion_types,
                mz_tolerance=mz_tolerance,
                add_losses=add_losses,
                add_first_prefix_ion=add_first_prefix_ion,
                remove_precursor_peaks=remove_precursor_peaks,
                normalize=normalize,
            )
            return df
        except Exception as e:
            error_msg = f"Annotation error ({spec_type} {spectrum.peptidoform}, charge={spectrum.precursor_charge}, peaks={len(spectrum.mz)}): {e}"
            results_summary["errors"].append(error_msg)
            add_log(f"    ❌ {error_msg}")
            return None

    def on_reference_loaded(variant_key: str, count: int):
        """Called when all references for a variant are loaded."""
        nonlocal refs_done, query_start_time
        add_log(f"[REF] {variant_key}: {count} reference spectra loaded")
        results_summary["reference_spectra"] += count

        # Update reference progress
        if expected_refs > 0:
            ref_frac = min(results_summary["reference_spectra"] / expected_refs, 1.0)
            ref_progress_bar.progress(ref_frac)

            elapsed = time.time() - ref_start_time
            if ref_frac > 0:
                estimated_total = elapsed / ref_frac
                remaining = estimated_total - elapsed
                ref_status.caption(f"Loaded {results_summary['reference_spectra']}/{expected_refs} references "
                                   f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)")

        # Check if all references are done
        if results_summary["reference_spectra"] >= expected_refs and not refs_done:
            refs_done = True
            if ref_progress_bar:
                ref_progress_bar.progress(1.0)
            ref_elapsed = time.time() - ref_start_time
            if ref_status:
                ref_status.caption(f"✅ Done: {results_summary['reference_spectra']} references loaded in {format_time(ref_elapsed)}")
            query_start_time = time.time()
            query_status.caption(f"Processing {expected_queries} query spectra...")

    def on_query_ready(result: StreamingResult):
        """Annotate query + build ComparisonTasks (no figure creation here)."""
        nonlocal refs_done, query_start_time
        query = result.query_spectrum
        refs = result.reference_spectra

        # Accumulate references (once per variant) and queries for cache
        if not use_parquet:
            if result.variant_key not in _cached_references:
                _cached_references[result.variant_key] = list(refs)
            _cached_queries.append((result.variant_key, query))
            n_queries = len(_cached_queries)
            if (n_queries - _last_incremental_save[0]) >= _INCREMENTAL_SAVE_INTERVAL:
                _last_incremental_save[0] = n_queries
                _try_save_cache(f"incremental {n_queries}", partial=True)

        # Mark refs as done if not already (first query means refs are done)
        if not refs_done:
            refs_done = True
            if ref_progress_bar:
                ref_progress_bar.progress(1.0)
            ref_elapsed = time.time() - ref_start_time
            if ref_status:
                ref_status.caption(f"✅ Done: {results_summary['reference_spectra']} references loaded in {format_time(ref_elapsed)}")
            query_start_time = time.time()

        progress_str = f"[{result.query_index + 1}/{result.total_queries}]"
        add_log(f"{progress_str} Annotating: {query.base_sequence} from {Path(query.source_file).name}")

        results_summary["query_spectra"] += 1

        # Update query progress with time estimate
        if expected_queries > 0:
            query_frac = min(results_summary["query_spectra"] / expected_queries, 1.0)
            query_progress_bar.progress(query_frac)

            if query_start_time:
                elapsed = time.time() - query_start_time
                if query_frac > 0:
                    estimated_total = elapsed / query_frac
                    remaining = estimated_total - elapsed
                    query_status.caption(f"Annotated {results_summary['query_spectra']}/{expected_queries} queries "
                                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)")

        # Annotate query spectrum
        query_df = annotate_spectrum(query, is_reference=False)
        if query_df is None:
            return

        # Ensure references are annotated (cache them)
        if result.variant_key not in annotated_references:
            annotated_references[result.variant_key] = []
            for ref in refs:
                ref_df = annotate_spectrum(ref, is_reference=True)
                if ref_df is not None:
                    annotated_references[result.variant_key].append((ref, ref_df))

        ref_data = annotated_references.get(result.variant_key, [])

        # Create folder for this variant
        if result.modification_key:
            folder_name = f"{result.base_sequence}_{sanitize_filename(result.modification_key)}"
            plot_title = f"{result.base_sequence} [{result.modification_key}]"
        else:
            folder_name = result.base_sequence
            plot_title = result.base_sequence

        if not ref_data:
            if not plot_unpaired:
                add_log(f"  Skipped (no reference): {Path(query.source_file).stem}")
                return

            # Build single spectrum task
            query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
            single_label = (f"{Path(query.source_file).stem}"
                            f"  (RT: {query_rt}, ID: {query.spectrum_id})")
            query_stem = shorten_stem(Path(query.source_file).stem)
            base_filename = f"{result.base_sequence}_{query_stem}"
            safe_filename = sanitize_filename(base_filename)

            comparison_tasks.append(ComparisonTask(
                task_id=task_counter[0],
                task_type="single",
                query_df=query_df,
                query_label=single_label,
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=str(png_dir),
                svg_dir=str(svg_dir),
                export_formats=export_formats,
                figsize=tuple(get_settings_value(_get_settings(), "export", "figsize", default=[12, 8])),
                dpi=get_settings_value(_get_settings(), "export", "dpi", default=300),
                normalize=normalize,
                query_file=Path(query.source_file).name,
                peptide=result.base_sequence,
                modification=result.modification_key or "",
                charge=query.precursor_charge,
            ))
            task_counter[0] += 1
            return

        # Build mirror task for each reference
        for ref, ref_df in ref_data:
            query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
            ref_rt = f"{ref.rt:.1f}" if ref.rt else "N/A"
            query_label = (f"Query: {Path(query.source_file).stem}"
                           f"  (RT: {query_rt}, ID: {query.spectrum_id})")
            ref_label = (f"Ref: {Path(ref.source_file).stem}"
                         f"  (RT: {ref_rt}, ID: {ref.spectrum_id})")

            query_stem = shorten_stem(Path(query.source_file).stem)
            ref_stem = shorten_stem(Path(ref.source_file).stem)
            base_filename = f"{result.base_sequence}_{query_stem}_vs_{ref_stem}"
            safe_filename = sanitize_filename(base_filename)

            comparison_tasks.append(ComparisonTask(
                task_id=task_counter[0],
                task_type="mirror",
                query_df=query_df,
                ref_df=ref_df,
                query_label=query_label,
                ref_label=ref_label,
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=str(png_dir),
                svg_dir=str(svg_dir),
                export_formats=export_formats,
                figsize=tuple(get_settings_value(_get_settings(), "export", "figsize", default=[12, 8])),
                dpi=get_settings_value(_get_settings(), "export", "dpi", default=300),
                normalize=normalize,
                query_file=Path(query.source_file).name,
                ref_file=Path(ref.source_file).name,
                peptide=result.base_sequence,
                modification=result.modification_key or "",
                charge=query.precursor_charge,
            ))
            task_counter[0] += 1

    def on_file_start(filename: str, is_reference: bool):
        """Called when starting to process a file."""
        file_type = "REFERENCE" if is_reference else "QUERY"
        add_log(f"Opening {file_type} file: {filename}")

    def on_file_complete(filename: str, loaded: int, not_found: int):
        """Called when done with a file."""
        add_log(f"Completed {filename}: {loaded} loaded, {not_found} not found")

    def on_error(error_msg: str):
        """Called on errors."""
        results_summary["errors"].append(error_msg)
        errors_metric.caption(f"⚠️ {len(results_summary['errors'])} errors")
        add_log(f"❌ {error_msg}")
    add_log(f"Export directory: {root_dir}")
    add_log(f"Formats: {export_formats}")
    add_log(f"Ion types: {ion_types}, m/z tolerance: {mz_tolerance} Da")
    add_log("-" * 50)

    try:
        streaming_workers = st.session_state.get(get_state_key('workers'), 4)

        if use_parquet:
            # FAST PATH: Load from Parquet
            add_log(f"Loading from Parquet: {parquet_file}")
            _parquet_refs = load_spectra_from_parquet(
                parquet_path=parquet_file,
                reference_substring=ref_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                error_callback=on_error,
            )
            _cached_references.update(_parquet_refs)
        elif partial_cache_path:
            # RESUME PATH: Load partial cache, then stream remaining from mzML
            add_log(f"Resuming from partial cache: {partial_cache_path}")

            # Step 1: Load partial Parquet through callbacks
            _parquet_refs = load_spectra_from_parquet(
                parquet_path=partial_cache_path,
                reference_substring=ref_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                error_callback=on_error,
            )
            _cached_references.update(_parquet_refs)

            cached_refs = results_summary["reference_spectra"]
            cached_queries = results_summary["query_spectra"]
            add_log(f"Partial cache loaded: {cached_refs} refs, {cached_queries} queries")

            # Step 2: Build skip set from partial cache
            from src.parquet_cache import get_cached_spectrum_ids
            skip_ids = get_cached_spectrum_ids(partial_cache_path)
            add_log(f"Skip set: {len(skip_ids)} spectra already cached")

            # Step 3: Stream remaining spectra from mzML
            _streaming_refs = load_spectra_streaming(
                variant_groups=variant_groups,
                mzml_paths=mzml_paths,
                reference_substring=ref_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                on_file_start=on_file_start,
                on_file_complete=on_file_complete,
                error_callback=on_error,
                workers=streaming_workers,
                skip_ids=skip_ids,
            )
            for vk, specs in _streaming_refs.items():
                if specs:
                    if vk not in _cached_references:
                        _cached_references[vk] = []
                    _cached_references[vk].extend(specs)
        else:
            # STANDARD PATH: Load from mzML
            _streaming_refs = load_spectra_streaming(
                variant_groups=variant_groups,
                mzml_paths=mzml_paths,
                reference_substring=ref_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                on_file_start=on_file_start,
                on_file_complete=on_file_complete,
                error_callback=on_error,
                workers=streaming_workers,
            )

        # Finalize annotation progress
        query_progress_bar.progress(1.0)
        annotation_time = time.time() - start_time
        query_status.caption(
            f"✅ Annotated {results_summary['query_spectra']} queries, "
            f"{len(comparison_tasks)} plots queued ({format_time(annotation_time)})"
        )

        # Save point: after streaming complete — final complete cache
        _try_save_cache("after loading", partial=False)

        # =================================================================
        # Parallel plot generation (with optional tmpfs staging)
        # =================================================================
        if comparison_tasks:
            from src.tmpfs_staging import TmpfsStaging, StagingConfig

            # Check for existing plot files (resumable runs)
            # Uses final disk paths — before any tmpfs rewriting
            skipped = 0
            for task in comparison_tasks:
                outputs_exist = True
                if "png" in task.export_formats and task.png_dir:
                    if not (Path(task.png_dir) / task.folder_name / f"{task.safe_filename}.png").exists():
                        outputs_exist = False
                if "svg" in task.export_formats and task.svg_dir:
                    if not (Path(task.svg_dir) / task.folder_name / f"{task.safe_filename}.svg").exists():
                        outputs_exist = False
                if outputs_exist:
                    task.skip_plot = True
                    skipped += 1

            if skipped > 0:
                add_log(f"Skipping {skipped}/{len(comparison_tasks)} plots (already on disk)")

            workers = st.session_state.get(get_state_key('workers'), 4)
            pending = len(comparison_tasks) - skipped
            st.markdown(f"**Generating plots** ({pending} new, {skipped} cached, {workers} workers)")
            plot_progress_bar = st.progress(0)
            plot_status = st.empty()
            plot_status.caption(f"Starting {len(comparison_tasks)} plot tasks...")

            plot_start = time.time()

            # Global counter for progress across chunks
            _global_completed = [0]
            total_tasks = len(comparison_tasks)

            def on_plot_complete_chunked(result, completed, total):
                _global_completed[0] += 1
                gc = _global_completed[0]
                if result.success:
                    results_summary["plots_generated"] += 1
                    if result.png_path:
                        results_summary["files_exported"]["png"].append(result.png_path)
                    if result.svg_path:
                        results_summary["files_exported"]["svg"].append(result.svg_path)
                    if result.summary_row:
                        summary_rows.append(result.summary_row)
                else:
                    results_summary["errors"].append(f"Plot error: {result.error}")

                frac = min(gc / total_tasks, 1.0)
                plot_progress_bar.progress(frac)
                elapsed = time.time() - plot_start
                if frac > 0:
                    remaining = (elapsed / frac) - elapsed
                    plot_status.caption(
                        f"Generated {gc}/{total_tasks} plots "
                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)"
                    )
                plots_metric.caption(f"📊 {results_summary['plots_generated']} plots generated")
                if results_summary["errors"]:
                    errors_metric.caption(f"⚠️ {len(results_summary['errors'])} errors")

            # Read tmpfs settings
            settings = _get_settings()
            staging_config = StagingConfig(
                enabled=get_settings_value(settings, "performance", "tmpfs_staging", default=True),
                tmpfs_base=get_settings_value(settings, "performance", "tmpfs_path", default="/dev/shm"),
                ram_fraction=get_settings_value(settings, "performance", "tmpfs_ram_fraction", default=0.40),
            )

            with TmpfsStaging(str(png_dir), str(svg_dir), export_formats,
                              staging_config, log_fn=add_log) as staging:
                for chunk in staging.iter_chunks(comparison_tasks):
                    run_parallel_comparisons(
                        tasks=chunk,
                        workers=workers,
                        on_task_complete=on_plot_complete_chunked,
                    )
                    staging.flush()

            plot_progress_bar.progress(1.0)
            plot_elapsed = time.time() - plot_start
            plot_status.caption(f"✅ Generated {results_summary['plots_generated']} plots in {format_time(plot_elapsed)}")

        # Save point 2: after plots
        _try_save_cache("after plots")

        # Final metrics
        total_time = time.time() - start_time
        plots_metric.caption(f"📊 {results_summary['plots_generated']} plots exported")
        if results_summary["errors"]:
            errors_metric.caption(f"⚠️ {len(results_summary['errors'])} errors")

        # Store results for display
        st.session_state[get_state_key('streaming_results')] = results_summary

        # Summary to terminal
        add_log("-" * 50)
        add_log("COMPLETE!")
        add_log(f"References: {results_summary['reference_spectra']}")
        add_log(f"Queries processed: {results_summary['query_spectra']}")
        add_log(f"Plots generated: {results_summary['plots_generated']}")
        add_log(f"Files exported: {len(results_summary['files_exported']['png'])} PNG, {len(results_summary['files_exported']['svg'])} SVG")
        add_log(f"Total time: {format_time(total_time)}")

        if results_summary["errors"]:
            add_log(f"Errors: {len(results_summary['errors'])}")

        st.success(f"✅ Complete! {results_summary['plots_generated']} plots exported to: {root_dir}")

        # =====================================================================
        # Write comparison summary CSV + XLSX
        # =====================================================================
        if summary_rows:
            try:
                import pandas as _pd
                summary_df = _pd.DataFrame(summary_rows)

                # Shorten filenames for readability
                all_filenames = list(set(
                    summary_df["query_file"].tolist() + summary_df["ref_file"].tolist()
                ))
                short_map = dict(zip(all_filenames, shorten_filenames(all_filenames)))
                summary_df.insert(0, "query_short", summary_df["query_file"].map(short_map))
                summary_df.insert(1, "ref_short", summary_df["ref_file"].map(short_map))

                # Reorder columns
                col_order = [
                    "query_short", "ref_short", "peptide", "modification", "charge",
                    "R2m", "R2b", "Cos", "Ent", "Dice", "Jacc", "ShR",
                    "total_q", "total_r", "ann_q", "ann_r", "shared",
                    "query_file", "ref_file",
                ]
                summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]

                # Rename for display
                rename_map = {
                    "query_short": "Query", "ref_short": "Reference",
                    "peptide": "Peptide", "modification": "Modification",
                    "charge": "Charge",
                    "R2m": "R² match", "R2b": "R² binned",
                    "Cos": "Cosine", "Ent": "Entropy",
                    "Dice": "Dice", "Jacc": "Jaccard", "ShR": "SharedRatio",
                    "total_q": "Peaks (Q)", "total_r": "Peaks (R)",
                    "ann_q": "Annotated (Q)", "ann_r": "Annotated (R)",
                    "shared": "Shared Ions",
                    "query_file": "Query File", "ref_file": "Ref File",
                }
                summary_df = summary_df.rename(columns=rename_map)

                # Write CSV (primary)
                csv_path = root_dir / "comparison_summary.csv"
                summary_df.to_csv(csv_path, index=False)
                results_summary["summary_csv"] = str(csv_path)
                add_log(f"Summary CSV: {csv_path}")

                # Write XLSX (secondary)
                try:
                    xlsx_path = root_dir / "comparison_summary.xlsx"
                    summary_df.to_excel(xlsx_path, index=False, sheet_name="Comparisons")
                    results_summary["summary_xlsx"] = str(xlsx_path)
                    add_log(f"Summary XLSX: {xlsx_path}")
                except Exception as xlsx_err:
                    add_log(f"Warning: Could not write XLSX: {xlsx_err}")

            except Exception as summary_err:
                add_log(f"Warning: Could not write summary: {summary_err}")
                logger.warning(f"Summary export failed: {summary_err}")

        # Save point 3: after summary export
        _try_save_cache("after summary")

    except Exception as e:
        st.error(f"Processing failed: {e}")
        logger.exception("Streaming processing error")


def display_streaming_results(results: dict):
    """Display streaming processing results."""
    refs = results.get("reference_spectra", 0)
    queries = results.get("query_spectra", 0)
    plots = results.get("plots_generated", 0)
    png_count = len(results.get("files_exported", {}).get("png", []))
    svg_count = len(results.get("files_exported", {}).get("svg", []))
    errors = results.get("errors", [])

    if refs > 0:
        st.caption(f"Processed {refs} reference + {queries} query spectra → {plots} mirror plots "
                   f"({png_count} PNG, {svg_count} SVG)")
    else:
        st.caption(f"Processed {queries} spectra → {plots} annotated spectrum plots "
                   f"({png_count} PNG, {svg_count} SVG)")

    # Show summary file paths
    summary_csv = results.get("summary_csv")
    summary_xlsx = results.get("summary_xlsx")
    if summary_csv or summary_xlsx:
        parts = []
        if summary_csv:
            parts.append(f"CSV: `{summary_csv}`")
        if summary_xlsx:
            parts.append(f"XLSX: `{summary_xlsx}`")
        st.caption(f"Comparison summary: {' | '.join(parts)}")

    if errors:
        with st.expander(f"⚠️ {len(errors)} Errors", expanded=False):
            for err in errors[:20]:
                st.caption(f"- {err}")
            if len(errors) > 20:
                st.caption(f"... and {len(errors) - 20} more")


def display_search_results(results: dict):
    """Display search results table."""
    from src.peptide_search import summarize_search_results

    summary_df = summarize_search_results(results)
    if not summary_df.empty:
        # Rename columns for display
        display_df = summary_df.rename(columns={
            'base_sequence': 'Peptide',
            'modification': 'Modification',
            'psm_count': 'PSMs',
            'idxml_files': 'idXML Files',
            'mzml_resolved': 'mzML Resolved',
            'unresolved': 'Unresolved'
        })
        # Drop display_label from table view
        display_df = display_df.drop(columns=['display_label'], errors='ignore')

        st.dataframe(
            display_df,
            width="stretch",
            hide_index=True,
            column_config={
                "Unresolved": st.column_config.NumberColumn(
                    "Unresolved",
                    help="PSMs that could not be mapped to mzML files"
                )
            }
        )


def render_results_gallery():
    """Render results summary."""
    st.header("Step 7: Results")

    streaming_results = st.session_state.get(get_state_key('streaming_results'))

    if not streaming_results:
        st.info("Run 'Process & Export' on the Analysis tab to generate results.")
        return

    st.success("✅ Processing Complete")

    display_streaming_results(streaming_results)

    export_root = streaming_results.get("export_root", "")
    if export_root:
        st.markdown("---")
        st.subheader("📁 Export Location")
        st.code(export_root)

        # Show exported files
        exported_files = streaming_results.get("files_exported", {})
        png_files = exported_files.get("png", [])
        svg_files = exported_files.get("svg", [])

        if png_files:
            with st.expander(f"PNG Files ({len(png_files)})", expanded=False):
                for f in png_files[:20]:
                    st.caption(f"📄 {Path(f).name}")
                if len(png_files) > 20:
                    st.caption(f"... and {len(png_files) - 20} more")

        if svg_files:
            with st.expander(f"SVG Files ({len(svg_files)})", expanded=False):
                for f in svg_files[:20]:
                    st.caption(f"📄 {Path(f).name}")
                if len(svg_files) > 20:
                    st.caption(f"... and {len(svg_files) - 20} more")


def render():
    """Main render function for the page."""
    init_session_state()

    st.title("Peptide Mirror Analysis")
    st.caption("Compare experimental spectra against reference spectra for target peptides")

    # Workflow steps as tabs
    steps = [
        ("1. Files", render_file_selection),
        ("2. Peptides", render_peptide_input),
        ("3. Reference", render_reference_selection),
        ("4. Settings", render_settings),
        ("5. Search", render_search),
        ("6. Process", render_analysis),
        ("7. Results", render_results_gallery),
    ]

    # Create tabs for steps
    tab_names = [s[0] for s in steps]
    tabs = st.tabs(tab_names)

    for tab, (name, render_func) in zip(tabs, steps):
        with tab:
            render_func()


# Run render when page is loaded
render()
