"""
Peptide Mirror Analysis Workspace

Interactive workspace for comparing query spectra against reference spectra.
"""
import streamlit as st
from pathlib import Path
from collections import defaultdict
import json
import logging

logger = logging.getLogger(__name__)

# State key prefix for this workspace
STATE_PREFIX = "peptide_mirror"

# ============================================================================
# SETTINGS LOADING - Shared with main.py
# ============================================================================

def load_settings() -> dict:
    """Load settings from settings.json (same logic as main.py)."""
    settings_path = Path(__file__).parent.parent / "settings.json"
    if settings_path.exists():
        with open(settings_path, 'r') as f:
            return json.load(f)
    return {}


def get_settings_value(settings: dict, *keys, default=None):
    """Safely get nested settings value."""
    value = settings
    for key in keys:
        if isinstance(value, dict):
            value = value.get(key)
        else:
            return default
    return value if value is not None else default


# Load settings once at module level
SETTINGS = load_settings()


def get_state_key(key: str) -> str:
    """Get namespaced session state key."""
    return f"{STATE_PREFIX}.{key}"


def init_session_state():
    """Initialize session state with defaults from settings.json."""
    # Get defaults from settings.json (same values used by main.py)
    ion_types_default = get_settings_value(SETTINGS, "ion_types", "default", default=["b", "y"])
    mz_tolerance_default = get_settings_value(SETTINGS, "mz_tolerance", "default", default=0.02)
    export_dir_default = get_settings_value(SETTINGS, "export", "default_directory", default="./exports")
    export_formats_default = get_settings_value(SETTINGS, "export", "formats", default=["png", "svg"])
    project_name_default = get_settings_value(SETTINGS, "export", "project_name", default="MSpecViz_analysis")
    data_dir_default = get_settings_value(SETTINGS, "data_directories", "default", default="/mnt/volume")

    # Annotation settings
    add_losses_default = get_settings_value(SETTINGS, "annotation", "add_losses", default=False)
    add_first_prefix_default = get_settings_value(SETTINGS, "annotation", "add_first_prefix_ion", default=True)
    remove_precursor_default = get_settings_value(SETTINGS, "annotation", "remove_precursor_peaks", default=False)

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

    # Show selection summary
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    if idxmls and mzmls:
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
            current_dir_selected_paths = [str(f) for f in filtered_files if f.name in selected_names]
            merged_selections = other_dir_selections + current_dir_selected_paths
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
                for folder, files in sorted(grouped.items()):
                    st.caption(f"**{folder}/** ({len(files)} files)")
                    for fname in sorted(files):
                        col1, col2 = st.columns([6, 1])
                        with col1:
                            st.text(f"  {fname}")
                        with col2:
                            full_path = str(Path(folder) / fname)
                            if st.button("❌", key=f"{key_prefix}_remove_{hash(full_path) % 100000}", help=f"Remove {fname}"):
                                updated = [p for p in all_selected if p != full_path]
                                st.session_state[get_state_key(state_key)] = updated
                                st.rerun()


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
    from core.peptide_search import get_unique_peptides_from_idxmls

    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    if not idxmls:
        st.warning("No idXML files selected")
        return

    # Get sample of peptides from files
    sample_peptides = get_unique_peptides_from_idxmls(idxmls, limit=50)

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
        "Reference identifier",
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
        st.info("Without a reference identifier, spectra will be compared query-to-query")
        return True

    return False


def render_settings():
    """Render analysis settings."""
    st.header("Step 4: Settings")

    # Get available ion types from settings
    available_ions = get_settings_value(SETTINGS, "ion_types", "available",
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

    # Info about precursor behavior
    if remove_precursor:
        st.caption("ℹ️ Precursor peaks will be removed from plots. Relative intensities will be recalculated.")
    elif "precursor" in st.session_state.get(get_state_key('ion_types'), []):
        st.caption("ℹ️ Precursor peaks will be shown and colored.")
    else:
        st.caption("ℹ️ Precursor peaks will be shown in grey (unannotated).")

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

    # Debug info
    with st.expander("Debug: Current State", expanded=False):
        st.write(f"idXML files: {len(idxmls)}")
        st.write(f"Peptide text length: {len(peptide_text)}")
        st.write(f"Peptides parsed: {peptides}")
        if idxmls:
            st.write(f"First idXML: {idxmls[0]}")

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
    else:
        st.caption("Search results will appear here after running the search.")


def render_analysis():
    """Render spectrum loading and analysis step."""
    st.header("Step 6: Process & Export")

    # Check prerequisites
    search_results = st.session_state.get(get_state_key('search_results'))
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])

    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')
    if not ref_substring:
        ref_substring = st.session_state.get(f"{STATE_PREFIX}_ref_input", '')

    if not search_results:
        st.warning("Please run peptide search first (Search tab)")
        return

    if not mzmls:
        st.warning("Please select mzML files first (Files tab)")
        return

    # Calculate expected spectra count
    total_psms = sum(len(vg.matches) for vg in search_results.values())
    streaming_results = st.session_state.get(get_state_key('streaming_results'))

    # Show status as caption (less prominent)
    if streaming_results:
        st.caption(f"✅ Processed {streaming_results.get('query_spectra', 0)} of {total_psms} PSMs, "
                   f"exported {streaming_results.get('plots_generated', 0)} plots")
    else:
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
                options=['png', 'svg'],
                default=st.session_state.get(get_state_key('export_formats'), ['png', 'svg']),
                key=f"{STATE_PREFIX}_export_formats_input"
            )
            st.session_state[get_state_key('export_formats')] = export_formats

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
    """Execute peptide search."""
    from core.peptide_search import search_peptides_in_idxmls

    progress_bar = st.progress(0, "Starting search...")
    status_text = st.empty()
    error_container = st.container()
    error_messages = []

    def progress_callback(current, total, message):
        if total > 0:
            progress_bar.progress(current / total, message)
            status_text.text(message)

    def error_callback(error_msg: str):
        """Handle errors immediately."""
        error_messages.append(error_msg)
        with error_container:
            st.warning(f"⚠️ {error_msg}")

    try:
        results = search_peptides_in_idxmls(
            idxmls, peptides,
            progress_callback=progress_callback,
            error_callback=error_callback
        )

        st.session_state[get_state_key('search_results')] = results
        progress_bar.progress(1.0, "Search complete!")
        status_text.text("Search complete!")

        if results:
            total_matches = sum(len(vg.matches) for vg in results.values())

            # Count resolved vs unresolved mzML mappings
            resolved = sum(1 for vg in results.values() for m in vg.matches if m.mzml_filename)
            unresolved = total_matches - resolved

            st.success(f"Found {total_matches} PSMs across {len(results)} peptide variants")

            if unresolved > 0:
                st.warning(f"{unresolved} PSMs could not be mapped to mzML files (missing id_merge_index)")

            # Show error summary
            if error_messages:
                with st.expander(f"⚠️ {len(error_messages)} warnings during search", expanded=False):
                    for msg in error_messages[:10]:
                        st.caption(msg)
                    if len(error_messages) > 10:
                        st.caption(f"... and {len(error_messages) - 10} more")
        else:
            st.warning("No matches found for the specified peptides")

    except Exception as e:
        st.error(f"Search failed: {e}")
        logger.exception("Peptide search error")


def run_streaming_processing():
    """
    Run streaming spectrum loading and processing.

    This mirrors the approach used in main.py:
    1. Load all reference spectra first
    2. Stream each query spectrum and immediately annotate, plot, and export
    """
    from core.spectrum_loader import load_spectra_streaming, StreamingResult
    from core.peak_annotation import prepare_spectrum_for_plot
    from core.mirror_plot import create_mirror_plot
    from core.export import sanitize_filename, create_export_directory
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Get state
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
    export_dir = st.session_state.get(get_state_key('export_dir'), './exports')
    export_formats = st.session_state.get(get_state_key('export_formats'), ['png', 'svg'])
    project_name = st.session_state.get(get_state_key('project_name'), 'MSpecViz_analysis')

    if not search_results:
        st.warning("No search results to process")
        return

    # Make a copy of search results
    import copy
    variant_groups = copy.deepcopy(search_results)

    mzml_paths = [Path(m) for m in mzmls]
    total_psms = sum(len(vg.matches) for vg in variant_groups.values())

    # Create export directory structure (same as main.py)
    root_dir, png_dir, svg_dir = create_export_directory(export_dir, project_name)

    # Count expected reference vs query PSMs
    ref_files_lower = {Path(m).name.lower() for m in mzmls
                       if ref_substring and ref_substring.lower() in Path(m).name.lower()}
    expected_refs = sum(
        1 for vg in variant_groups.values() for m in vg.matches
        if m.mzml_filename and m.mzml_filename.lower() in ref_files_lower
    )
    expected_queries = total_psms - expected_refs

    # Helper to format time as hours/minutes/seconds
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

    # Progress display - two separate progress bars
    import time
    start_time = time.time()
    ref_start_time = start_time
    query_start_time = None

    # Reference progress section
    st.markdown("**Reference Spectra**")
    ref_progress_bar = st.progress(0)
    ref_status = st.empty()
    ref_status.caption(f"Loading {expected_refs} reference spectra...")

    # Query progress section
    st.markdown("**Query Spectra**")
    query_progress_bar = st.progress(0)
    query_status = st.empty()
    query_status.caption(f"Waiting for references to complete...")

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
    refs_done = False

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
            ref_progress_bar.progress(1.0)
            ref_elapsed = time.time() - ref_start_time
            ref_status.caption(f"✅ Done: {results_summary['reference_spectra']} references loaded in {format_time(ref_elapsed)}")
            query_start_time = time.time()
            query_status.caption(f"Processing {expected_queries} query spectra...")

    def on_query_ready(result: StreamingResult):
        """Called for each query spectrum - annotate, plot, export immediately (same as main.py)."""
        nonlocal refs_done, query_start_time
        query = result.query_spectrum
        refs = result.reference_spectra

        # Mark refs as done if not already (first query means refs are done)
        if not refs_done:
            refs_done = True
            ref_progress_bar.progress(1.0)
            ref_elapsed = time.time() - ref_start_time
            ref_status.caption(f"✅ Done: {results_summary['reference_spectra']} references loaded in {format_time(ref_elapsed)}")
            query_start_time = time.time()

        progress_str = f"[{result.query_index + 1}/{result.total_queries}]"
        add_log(f"{progress_str} Processing: {query.base_sequence} from {Path(query.source_file).name}")

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
                    query_status.caption(f"Processed {results_summary['query_spectra']}/{expected_queries} queries "
                                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)")

        # Annotate query spectrum
        query_df = annotate_spectrum(query, is_reference=False)
        if query_df is None:
            # Error already logged in annotate_spectrum
            return

        # Ensure references are annotated (cache them)
        if result.variant_key not in annotated_references:
            annotated_references[result.variant_key] = []
            for ref in refs:
                ref_df = annotate_spectrum(ref, is_reference=True)
                if ref_df is not None:
                    annotated_references[result.variant_key].append((ref, ref_df))

        ref_data = annotated_references.get(result.variant_key, [])

        if not ref_data:
            add_log(f"  ⚠️ No annotated references available")
            return

        # Create folder for this variant
        if result.modification_key:
            folder_name = f"{result.base_sequence}_{sanitize_filename(result.modification_key)}"
        else:
            folder_name = result.base_sequence

        # Generate mirror plot for each reference (same as main.py)
        for ref, ref_df in ref_data:
            try:
                # Create plot
                fig = create_mirror_plot(
                    upper_spectrum=query_df,
                    lower_spectrum=ref_df,
                    upper_label=f"Query: {Path(query.source_file).stem}",
                    lower_label=f"Ref: {Path(ref.source_file).stem}",
                    title=f"{result.base_sequence}",
                )

                results_summary["plots_generated"] += 1
                plots_metric.caption(f"📊 {results_summary['plots_generated']} plots generated")

                # Generate filename (same as main.py)
                query_stem = Path(query.source_file).stem[:30]
                ref_stem = Path(ref.source_file).stem[:30]
                base_filename = f"{result.base_sequence}_{query_stem}_vs_{ref_stem}"
                safe_filename = sanitize_filename(base_filename)

                # Export PNG
                if "png" in export_formats:
                    variant_png_dir = png_dir / folder_name
                    variant_png_dir.mkdir(exist_ok=True)
                    png_path = variant_png_dir / f"{safe_filename}.png"
                    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
                    results_summary["files_exported"]["png"].append(str(png_path))

                # Export SVG
                if "svg" in export_formats:
                    variant_svg_dir = svg_dir / folder_name
                    variant_svg_dir.mkdir(exist_ok=True)
                    svg_path = variant_svg_dir / f"{safe_filename}.svg"
                    fig.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
                    results_summary["files_exported"]["svg"].append(str(svg_path))

                plt.close(fig)
                add_log(f"  ✓ Exported: {safe_filename}")

            except Exception as e:
                results_summary["errors"].append(f"Plot/export error: {e}")
                errors_metric.caption(f"⚠️ {len(results_summary['errors'])} errors")
                add_log(f"  ❌ Error: {e}")

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
        # Run the streaming pipeline (same as main.py)
        reference_spectra = load_spectra_streaming(
            variant_groups=variant_groups,
            mzml_paths=mzml_paths,
            reference_substring=ref_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=on_file_complete,
            error_callback=on_error,
        )

        # Finalize progress bars
        query_progress_bar.progress(1.0)
        total_time = time.time() - start_time
        query_status.caption(f"✅ Done: {results_summary['query_spectra']} queries processed in {format_time(total_time)} total")

        # Final metrics
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

    except Exception as e:
        st.error(f"Processing failed: {e}")
        status_text.text(f"Error: {e}")
        logger.exception("Streaming processing error")


def display_streaming_results(results: dict):
    """Display streaming processing results."""
    refs = results.get("reference_spectra", 0)
    queries = results.get("query_spectra", 0)
    plots = results.get("plots_generated", 0)
    png_count = len(results.get("files_exported", {}).get("png", []))
    svg_count = len(results.get("files_exported", {}).get("svg", []))
    errors = results.get("errors", [])

    st.caption(f"Processed {refs} reference + {queries} query spectra → {plots} plots "
               f"({png_count} PNG, {svg_count} SVG)")

    if errors:
        with st.expander(f"⚠️ {len(errors)} Errors", expanded=False):
            for err in errors[:20]:
                st.caption(f"- {err}")
            if len(errors) > 20:
                st.caption(f"... and {len(errors) - 20} more")


def display_search_results(results: dict):
    """Display search results table."""
    from core.peptide_search import summarize_search_results

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
    """Main render function for the workspace."""
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
