"""
Peptide Mirror Analysis

Interactive page for exploring peptide-spectrum matches.
Combines openms_insight interactive components (Table, LinePlot, SequenceView,
Heatmap) with Plotly mirror plots, full MSpecViz comparison metrics,
Koina AI predictions, parallel export, and export history.
"""
import io
import json
import os
import shutil
import zipfile
import streamlit as st
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import logging
import polars as pl

from src.common.common import page_setup
from utils.helpers import format_time, get_nested as get_settings_value, parse_peptide_list

params = page_setup()

logger = logging.getLogger(__name__)

# State key prefix for this workspace
STATE_PREFIX = "peptide_mirror"


def _get_settings() -> dict:
    """Get settings from session state (loaded by page_setup)."""
    return st.session_state.get("settings", {})


def get_state_key(key: str) -> str:
    """Get namespaced session state key."""
    return f"{STATE_PREFIX}.{key}"


def init_session_state(params: dict = None):
    """Initialize session state with defaults from settings.json."""
    if params:
        for key, value in params.get("peptide_mirror", {}).items():
            if key.startswith("peptide_mirror.") and key not in st.session_state:
                st.session_state[key] = value

    _settings = _get_settings()
    ion_types_default = get_settings_value(_settings, "ion_types", "default", default=["b", "y"])
    mz_tolerance_default = get_settings_value(_settings, "mz_tolerance", "default", default=0.02)
    workspace = st.session_state.get('workspace', None)
    export_dir_default = (
        str(Path(workspace) / "exports")
        if workspace
        else get_settings_value(_settings, "export", "default_directory", default="./exports")
    )
    add_losses_default = get_settings_value(_settings, "annotation", "add_losses", default=False)
    add_first_prefix_default = get_settings_value(_settings, "annotation", "add_first_prefix_ion", default=True)
    remove_precursor_default = get_settings_value(_settings, "annotation", "remove_precursor_peaks", default=False)
    workers_default = get_settings_value(_settings, "performance", "workers", default=4)

    defaults = {
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
        'workers': workers_default,
        'intensity_mode': 'Relative',
        'variant_groups': None,
        'data_ready': False,
        'koina_enabled': False,
    }

    for key, default in defaults.items():
        state_key = get_state_key(key)
        if state_key not in st.session_state:
            st.session_state[state_key] = default


# ---------------------------------------------------------------------------
# Cached filesystem helpers
# ---------------------------------------------------------------------------

@st.cache_data(ttl=10)
def _list_subdirectories(path: str) -> list:
    try:
        return sorted(str(d) for d in Path(path).iterdir() if d.is_dir())
    except (PermissionError, OSError):
        return []


@st.cache_data(ttl=10)
def _list_files(path: str, extensions: tuple, recursive: bool = False, _max_files: int = 5000) -> list:
    try:
        files = []
        for ext in extensions:
            gen = Path(path).rglob(f"*{ext}") if recursive else Path(path).glob(f"*{ext}")
            for f in gen:
                files.append(f)
                if len(files) >= _max_files:
                    break
            if len(files) >= _max_files:
                break
        return sorted(str(f) for f in files)
    except (PermissionError, OSError):
        return []


# ---------------------------------------------------------------------------
# Tab 1: Setup
# ---------------------------------------------------------------------------

def render_directory_browser(key_prefix: str, extensions: list, state_key: str, default_dir: str = "/"):
    """Render a directory browser with file selection."""
    dir_key = f"{key_prefix}_dir"
    nav_key = f"{key_prefix}_nav_to"
    select_all_key = f"{key_prefix}_select_all"
    recursive_key = f"{key_prefix}_recursive"

    if dir_key not in st.session_state:
        st.session_state[dir_key] = default_dir

    if nav_key in st.session_state and st.session_state[nav_key]:
        st.session_state[dir_key] = st.session_state[nav_key]
        st.session_state[nav_key] = None
        if select_all_key in st.session_state:
            del st.session_state[select_all_key]
        if recursive_key in st.session_state:
            del st.session_state[recursive_key]

    current_dir = st.session_state[dir_key]

    col1, col2 = st.columns([6, 1])
    with col1:
        new_path = st.text_input(
            "Path", value=current_dir,
            key=f"{key_prefix}_path_input_{abs(hash(current_dir))}",
            label_visibility="collapsed", placeholder="Enter path..."
        )
        if new_path != current_dir and Path(new_path).exists() and Path(new_path).is_dir():
            st.session_state[dir_key] = new_path
            if select_all_key in st.session_state:
                del st.session_state[select_all_key]
            if recursive_key in st.session_state:
                del st.session_state[recursive_key]
            st.rerun()

    with col2:
        parent = str(Path(current_dir).parent)
        if st.button("up", key=f"{key_prefix}_up_btn", help="Parent directory",
                      disabled=(parent == current_dir)):
            st.session_state[nav_key] = parent
            st.rerun()

    try:
        subdir_paths = _list_subdirectories(current_dir)
        subdirs = [Path(d) for d in subdir_paths]
        if subdirs:
            st.markdown(
                "<style>.folder-nav button{font-size:0.8rem !important;"
                "padding:0.15rem 0.4rem !important;min-height:0 !important}</style>",
                unsafe_allow_html=True,
            )
            cols_per_row = 4
            max_folders = 16
            for i in range(0, min(len(subdirs), max_folders), cols_per_row):
                cols = st.columns(cols_per_row)
                for j, col in enumerate(cols):
                    idx = i + j
                    if idx < len(subdirs) and idx < max_folders:
                        folder = subdirs[idx]
                        with col:
                            if st.button(f"{folder.name}", key=f"{key_prefix}_folder_{idx}",
                                         help=str(folder)):
                                st.session_state[nav_key] = str(folder)
                                st.rerun()
            if len(subdirs) > max_folders:
                st.caption(f"... and {len(subdirs) - max_folders} more folders")
    except Exception as e:
        st.error(f"Cannot access: {e}")

    recursive = st.checkbox("Include subdirectories", value=False,
                            key=f"{key_prefix}_recursive",
                            help="Scan all subdirectories for matching files")

    try:
        all_files = [Path(f) for f in _list_files(current_dir, tuple(extensions), recursive)]
        current_path = Path(current_dir)
        file_labels = ([str(f.relative_to(current_path)) for f in all_files]
                       if recursive else [f.name for f in all_files])
        label_to_path = {lbl: str(f) for lbl, f in zip(file_labels, all_files)}

        filter_text = st.text_input("Filter files", key=f"{key_prefix}_filter",
                                    placeholder="Type to filter filenames...")
        if filter_text:
            filter_lower = filter_text.lower()
            filtered_labels = [lbl for lbl in file_labels if filter_lower in lbl.lower()]
        else:
            filtered_labels = file_labels

        if filtered_labels:
            all_selected_paths = set(st.session_state.get(get_state_key(state_key), []))
            current_scope_paths = {label_to_path[lbl] for lbl in file_labels}
            if recursive:
                other_dir_selections = [p for p in st.session_state.get(get_state_key(state_key), [])
                                        if p not in current_scope_paths]
            else:
                other_dir_selections = [p for p in st.session_state.get(get_state_key(state_key), [])
                                        if Path(p).parent != current_path]

            valid_defaults = [lbl for lbl in filtered_labels if label_to_path[lbl] in all_selected_paths]
            all_in_scope_selected = set(valid_defaults) == set(filtered_labels) and len(valid_defaults) > 0

            select_all = st.checkbox(
                f"Select all {len(filtered_labels)} files" +
                (f" (filtered from {len(file_labels)})" if filter_text else ""),
                value=all_in_scope_selected, key=select_all_key
            )

            if select_all:
                selected_labels = filtered_labels
            else:
                selected_labels = st.multiselect(
                    "Select files" + (" (recursive):" if recursive else ":"),
                    options=filtered_labels,
                    default=valid_defaults if not all_in_scope_selected else [],
                    key=f"{key_prefix}_files"
                )

            current_selected_paths = [label_to_path[lbl] for lbl in selected_labels]
            merged_selections = list(dict.fromkeys(other_dir_selections + current_selected_paths))
            st.session_state[get_state_key(state_key)] = merged_selections

        else:
            if filter_text:
                st.info(f"No files match filter '{filter_text}'")
            else:
                st.info(f"No {extensions[0]} files found")

    except Exception as e:
        st.error(f"Error listing files: {e}")

    # Clear selections button
    all_selected = st.session_state.get(get_state_key(state_key), [])
    if all_selected:
        if st.button(f"Clear all selections ({len(all_selected)} files)",
                     key=f"{key_prefix}_clear_all", width="stretch"):
            st.session_state[get_state_key(state_key)] = []
            if select_all_key in st.session_state:
                del st.session_state[select_all_key]
            st.rerun()


def render_peptide_input():
    st.subheader("Target Peptides")
    col1, col2 = st.columns([2, 1])

    with col1:
        peptide_text = st.text_area(
            "Enter peptide sequences (one per line)",
            value=st.session_state.get(get_state_key('target_peptides'), ''),
            height=150, key=f"{STATE_PREFIX}_peptide_input",
            placeholder="PEPTIDESEQ\nANOTHERPEPTIDE\n..."
        )
        st.session_state[get_state_key('target_peptides')] = peptide_text

    with col2:
        st.markdown("**Tips:**\n- One peptide per line\n"
                     "- Without modifications\n"
                     "- Case insensitive")
        peptides = parse_peptide_list(peptide_text)
        if peptides:
            st.success(f"{len(peptides)} peptide(s) entered")
    return len(peptides) > 0


def render_reference_selection():
    st.subheader("Reference Identification")
    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    if not mzmls:
        st.warning("Please select mzML files first")
        return False

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")
    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=st.session_state.get(get_state_key('reference_substring'), ''),
        key=f"{STATE_PREFIX}_ref_input", placeholder="e.g., Synpep"
    )
    st.session_state[get_state_key('reference_substring')] = ref_substring

    if ref_substring:
        reference_files = [Path(m).name for m in mzmls if ref_substring.lower() in Path(m).name.lower()]
        query_files = [Path(m).name for m in mzmls if ref_substring.lower() not in Path(m).name.lower()]

        with st.expander(f"Reference ({len(reference_files)}) / Query ({len(query_files)}) files", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Reference files ({len(reference_files)})**")
                for f in reference_files[:5]:
                    st.caption(f)
            with col2:
                st.markdown(f"**Query files ({len(query_files)})**")
                for f in query_files[:5]:
                    st.caption(f)

        if reference_files and query_files:
            st.success("Reference/Query separation configured")
            return True
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")
    else:
        st.info("No reference identifier -- individual annotated spectrum plots only.")
        return True

    return False


def render_settings():
    st.subheader("Settings")
    available_ions = get_settings_value(_get_settings(), "ion_types", "available",
                                        default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Ion Types**")
        current_ions = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])
        selected_ions = st.multiselect("Select ion types", options=available_ions,
                                       default=[i for i in current_ions if i in available_ions],
                                       key=f"{STATE_PREFIX}_ion_select")
        st.session_state[get_state_key('ion_types')] = selected_ions

    with col2:
        st.markdown("**m/z Tolerance**")
        tolerance = st.number_input("Tolerance (Da)", min_value=0.001, max_value=1.0,
                                    value=float(st.session_state.get(get_state_key('mz_tolerance'), 0.02)),
                                    step=0.005, format="%.3f", key=f"{STATE_PREFIX}_tol_input")
        st.session_state[get_state_key('mz_tolerance')] = tolerance

    st.markdown("**Annotation Settings**")
    ann_col1, ann_col2, ann_col3 = st.columns(3)

    with ann_col1:
        add_losses = st.checkbox("Add neutral losses (H2O, NH3)",
                                 value=st.session_state.get(get_state_key('add_losses'), False),
                                 key=f"{STATE_PREFIX}_add_losses")
        st.session_state[get_state_key('add_losses')] = add_losses

    with ann_col2:
        add_first_prefix = st.checkbox("Add first prefix ion",
                                       value=st.session_state.get(get_state_key('add_first_prefix_ion'), True),
                                       key=f"{STATE_PREFIX}_add_first_prefix")
        st.session_state[get_state_key('add_first_prefix_ion')] = add_first_prefix

    with ann_col3:
        remove_precursor = st.checkbox("Remove precursor peaks",
                                       value=st.session_state.get(get_state_key('remove_precursor_peaks'), False),
                                       key=f"{STATE_PREFIX}_remove_precursor")
        st.session_state[get_state_key('remove_precursor_peaks')] = remove_precursor

    intensity_mode = st.radio("**Intensity Display Mode**",
                              options=["Relative", "Raw"],
                              captions=["Relative to the highest intensity", "Raw intensity"],
                              index=0 if st.session_state.get(get_state_key('intensity_mode'), 'Relative') == 'Relative' else 1,
                              horizontal=True, key=f"{STATE_PREFIX}_intensity_mode",
                              label_visibility="visible")
    st.session_state[get_state_key('intensity_mode')] = intensity_mode

    st.markdown("**Performance**")
    workers = st.slider("Parallel workers", min_value=1, max_value=32,
                        value=st.session_state.get(get_state_key('workers'), 4),
                        key=f"{STATE_PREFIX}_workers",
                        help="Number of parallel processes for plot generation (~30 MB per worker)")
    st.session_state[get_state_key('workers')] = workers

    # --- Koina AI Prediction Settings ---
    st.markdown("---")
    st.markdown("**Koina AI Prediction**")
    st.caption(
        "When enabled, spectra without an experimental reference are compared "
        "against AI-predicted fragmentation patterns (intensity + m/z)."
    )

    try:
        from src.koina_prediction import (
            KOINA_AVAILABLE, KOINA_MODELS, INSTRUMENT_TYPES,
            DEFAULT_MODEL, DEFAULT_COLLISION_ENERGY, DEFAULT_INSTRUMENT_TYPE,
            get_koina_settings,
        )
    except ImportError:
        KOINA_AVAILABLE = False
        KOINA_MODELS = {}
        INSTRUMENT_TYPES = []
        DEFAULT_MODEL = ""
        DEFAULT_COLLISION_ENERGY = 35
        DEFAULT_INSTRUMENT_TYPE = "TIMSTOF"

        def get_koina_settings(s):
            return {"enabled": False, "model": DEFAULT_MODEL,
                    "collision_energy": DEFAULT_COLLISION_ENERGY,
                    "instrument_type": DEFAULT_INSTRUMENT_TYPE}

    _k_defaults = get_koina_settings(_get_settings())

    koina_enabled = st.checkbox(
        "Enable Koina fallback reference",
        value=st.session_state.get(get_state_key("koina_enabled"), _k_defaults["enabled"]),
        key=f"{STATE_PREFIX}_koina_enabled_cb",
        help="Use AI-predicted spectrum as reference when no experimental reference mzML exists",
        disabled=not KOINA_AVAILABLE,
    )
    st.session_state[get_state_key("koina_enabled")] = koina_enabled

    if not KOINA_AVAILABLE:
        st.warning("koinapy is not installed. Install it to enable Koina predictions.")

    if koina_enabled and KOINA_AVAILABLE:
        k_col1, k_col2, k_col3 = st.columns(3)

        with k_col1:
            model_names = list(KOINA_MODELS.keys())
            model_ids = list(KOINA_MODELS.values())
            current_model_id = st.session_state.get(
                get_state_key("koina_model"), _k_defaults["model"],
            )
            try:
                default_idx = model_ids.index(current_model_id)
            except ValueError:
                default_idx = 0

            selected_model_name = st.selectbox(
                "Prediction Model", options=model_names, index=default_idx,
                key=f"{STATE_PREFIX}_koina_model_select",
                help="Deep-learning model for fragment ion prediction",
            )
            st.session_state[get_state_key("koina_model")] = KOINA_MODELS[selected_model_name]

        with k_col2:
            ce = st.number_input(
                "Collision Energy (eV)", min_value=1, max_value=100,
                value=st.session_state.get(get_state_key("koina_ce"), _k_defaults["collision_energy"]),
                key=f"{STATE_PREFIX}_koina_ce_input",
            )
            st.session_state[get_state_key("koina_ce")] = ce

        with k_col3:
            current_inst = st.session_state.get(
                get_state_key("koina_instrument"), _k_defaults["instrument_type"],
            )
            try:
                inst_idx = INSTRUMENT_TYPES.index(current_inst)
            except ValueError:
                inst_idx = 0

            instrument = st.selectbox(
                "Instrument Type", options=INSTRUMENT_TYPES, index=inst_idx,
                key=f"{STATE_PREFIX}_koina_inst_select",
            )
            st.session_state[get_state_key("koina_instrument")] = instrument

    return True


@st.fragment
def render_setup_tab():
    _settings = _get_settings()
    data_dir_default = get_settings_value(_settings, "data_directories", "default", default="/")

    st.header("File Selection")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_idxml",
                                 extensions=['.idXML', '.idxml'],
                                 state_key='selected_idxmls',
                                 default_dir=data_dir_default)
    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_mzml",
                                 extensions=['.mzML', '.mzml'],
                                 state_key='selected_mzmls',
                                 default_dir=data_dir_default)

    st.markdown("---")
    render_peptide_input()
    st.markdown("---")
    render_reference_selection()
    st.markdown("---")
    render_settings()


@st.cache_data(show_spinner="Scanning idXML files...")
def _load_idxml_peptides(idxml_paths: tuple) -> list:
    """Load unique peptide sequences from idXML files (cached)."""
    import pyopenms as oms
    all_peptides = set()
    for idxml_path in idxml_paths[:5]:
        prot_ids = []
        pep_ids = []
        oms.IdXMLFile().load(idxml_path, prot_ids, pep_ids)
        for pid in pep_ids:
            for hit in pid.getHits():
                seq = str(hit.getSequence().toUnmodifiedString())
                if seq:
                    all_peptides.add(seq)
    return sorted(all_peptides)


def _render_idxml_peptide_preview(idxmls: list):
    """Show unique peptide sequences found in selected idXML files."""
    try:
        sorted_peps = _load_idxml_peptides(tuple(idxmls))
        if sorted_peps:
            st.caption(f"{len(sorted_peps)} unique peptides found in {min(len(idxmls), 5)} file(s)")
            st.text_area("Peptides found", value="\n".join(sorted_peps[:100]),
                         height=150, disabled=True, key=f"{STATE_PREFIX}_pep_preview")
            if len(sorted_peps) > 100:
                st.caption(f"... and {len(sorted_peps) - 100} more")
        else:
            st.info("No peptide hits found in selected files")
    except Exception as e:
        st.warning(f"Could not preview peptides: {e}")


# ---------------------------------------------------------------------------
# Tab 2: Search & Load
# ---------------------------------------------------------------------------

def run_peptide_search(idxmls: list, peptides: list):
    """Execute peptide search."""
    from src.peptide_search import search_peptides_in_idxmls
    import time as _time

    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    error_messages = []

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
                    status_text.caption(f"{current}/{total} files  "
                                        f"({format_time(elapsed)} elapsed, ~{format_time(remaining)} remaining)")
                elif current >= total:
                    status_text.caption(f"Done in {format_time(elapsed)}")

        def error_callback(error_msg: str):
            error_messages.append(error_msg)

        try:
            workers = st.session_state.get(get_state_key('workers'), 4)
            results = search_peptides_in_idxmls(
                idxmls, peptides, mzml_paths=mzmls,
                progress_callback=progress_callback, error_callback=error_callback,
                workers=workers,
            )
            st.session_state[get_state_key('variant_groups')] = results
            progress_bar.progress(1.0, "Search complete!")

            if results:
                total_matches = sum(len(vg.matches) for vg in results.values())
                status.update(label=f"Found {total_matches} PSMs across {len(results)} variants",
                              state="complete")
                from src.common.common import save_params
                save_params(params)

                unresolved = sum(1 for vg in results.values() for m in vg.matches if not m.mzml_filename)
                if unresolved > 0:
                    st.warning(f"{unresolved} PSMs could not be mapped to mzML files")
                if error_messages:
                    with st.expander(f"{len(error_messages)} warnings", expanded=False):
                        for msg in error_messages[:10]:
                            st.caption(msg)
            else:
                status.update(label="No matches found", state="error")
        except Exception as e:
            status.update(label=f"Search failed: {e}", state="error")
            logger.exception("Peptide search error")


def run_search_and_load():
    """Stream spectra, capture references, annotate, compute metrics, build DataFrames."""
    from pathlib import Path as _Path
    import copy
    import time as _time
    from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key
    from src.peak_annotation import prepare_spectrum_for_plot
    from src.parquet_builder import (
        build_psm_table, build_peaks_dataframe,
        build_sequence_dataframe, build_heatmap_dataframe,
    )
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text, METRIC_SHORT_TO_DISPLAY

    variant_groups = st.session_state.get(get_state_key('variant_groups'))
    if not variant_groups:
        st.error("No search results. Run peptide search first.")
        return

    mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
    if not mzmls:
        st.error("No mzML files selected.")
        return

    ref_substring = st.session_state.get(get_state_key('reference_substring'), '')
    ion_types = st.session_state.get(get_state_key('ion_types'), ['b', 'y'])
    mz_tolerance = st.session_state.get(get_state_key('mz_tolerance'), 0.02)
    normalize = st.session_state.get(get_state_key('intensity_mode'), 'Relative') == 'Relative'
    workers = st.session_state.get(get_state_key('workers'), 4)

    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)

    # Storage
    annotated_queries = []      # (psm_id, LoadedSpectrum, annotated_df)
    reference_store = {}        # variant_key -> [(LoadedSpectrum, annotated_df)]
    comparison_rows = []        # dicts for comparison_summary

    total_psms = sum(len(vg.matches) for vg in vg_copy.values())

    add_losses = st.session_state.get(get_state_key('add_losses'), False)
    add_first_prefix_ion = st.session_state.get(get_state_key('add_first_prefix_ion'), True)
    remove_precursor_peaks = st.session_state.get(get_state_key('remove_precursor_peaks'), False)

    annotation_kwargs = dict(
        ion_types=ion_types,
        mz_tolerance=mz_tolerance,
        normalize=normalize,
        add_losses=add_losses,
        add_first_prefix_ion=add_first_prefix_ion,
        remove_precursor_peaks=remove_precursor_peaks,
    )

    with st.status(f"Loading {total_psms} spectra...", expanded=True) as status:
        progress_bar = st.progress(0)
        status_text = st.empty()
        _load_start = _time.time()
        loaded_count = [0]

        def on_reference_loaded(variant_key, count):
            pass

        def on_query_ready(result):
            """Process each query: annotate, pair with references, compute metrics."""
            key = spectrum_lookup_key(result.query_spectrum)
            psm_id = psm_lookup.get(key)
            if psm_id is None:
                return

            try:
                query_df = prepare_spectrum_for_plot(
                    result.query_spectrum.mz,
                    result.query_spectrum.intensity,
                    result.query_spectrum.peptidoform,
                    result.query_spectrum.precursor_charge,
                    **annotation_kwargs,
                )
                annotated_queries.append((psm_id, result.query_spectrum, query_df))
            except Exception as e:
                logger.warning(f"Annotation failed psm_id={psm_id}: {e}")
                return

            # Annotate and cache references for this variant (once)
            vk = result.variant_key
            if vk not in reference_store and result.reference_spectra:
                ref_list = []
                for ref_spec in result.reference_spectra:
                    try:
                        ref_df = prepare_spectrum_for_plot(
                            ref_spec.mz, ref_spec.intensity,
                            ref_spec.peptidoform, ref_spec.precursor_charge,
                            **annotation_kwargs,
                        )
                        ref_list.append((ref_spec, ref_df))
                    except Exception as e:
                        logger.warning(f"Reference annotation failed: {e}")
                reference_store[vk] = ref_list

            # Koina fallback: if no experimental references, try AI prediction
            refs = reference_store.get(vk, [])
            if not refs and st.session_state.get(get_state_key('koina_enabled'), False):
                koina_ref = _get_koina_prediction(
                    result.query_spectrum.peptidoform,
                    result.query_spectrum.precursor_charge,
                    annotation_kwargs,
                )
                if koina_ref is not None:
                    reference_store[vk] = [koina_ref]
                    refs = reference_store[vk]

            # Compute comparison metrics against each reference
            for ref_spec, ref_df in refs:
                try:
                    metrics = compute_all_metrics(query_df, ref_df)
                    row = {
                        "psm_id": psm_id,
                        "Peptide": result.base_sequence,
                        "Modification": result.modification_key,
                        "Charge": result.query_spectrum.precursor_charge,
                        "Query": _Path(result.query_spectrum.source_file).stem,
                        "Reference": _Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina",
                        "query_spectrum_id": result.query_spectrum.spectrum_id,
                        "ref_spectrum_id": ref_spec.spectrum_id,
                        "query_rt": result.query_spectrum.rt,
                        "ref_rt": ref_spec.rt,
                    }
                    for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
                        row[display_name] = metrics.get(short_key, float("nan"))
                    for count_key in ("total_q", "total_r", "ann_q", "ann_r", "shared"):
                        row[count_key] = metrics.get(count_key, 0)

                    # Ion match statistics (17 supplementary stats)
                    try:
                        from src.parallel_processor import compute_ion_match_statistics
                        ims_stats = compute_ion_match_statistics(
                            query_df, result.query_spectrum.peptidoform,
                            result.query_spectrum.precursor_charge,
                            result.query_spectrum.precursor_mz,
                        )
                        row.update(ims_stats)
                    except Exception:
                        pass

                    comparison_rows.append(row)
                except Exception as e:
                    logger.warning(f"Metric computation failed: {e}")

            loaded_count[0] += 1
            if loaded_count[0] % 10 == 0 or loaded_count[0] == total_psms:
                frac = min(loaded_count[0] / max(total_psms, 1), 1.0)
                elapsed = _time.time() - _load_start
                progress_bar.progress(frac)
                status_text.caption(f"{loaded_count[0]}/{total_psms} spectra "
                                    f"({format_time(elapsed)} elapsed)")

        def on_file_start(path, is_reference):
            kind = "reference" if is_reference else "query"
            status_text.caption(f"Processing {kind}: {_Path(path).name}")

        def on_file_complete(path, loaded, total):
            logger.info(f"Completed {_Path(path).name}: {loaded}/{total} spectra")

        def error_callback(error):
            logger.error(f"Streaming error: {error}")

        load_spectra_streaming(
            variant_groups=vg_copy,
            mzml_paths=[_Path(m) for m in mzmls],
            reference_substring=ref_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=on_file_complete,
            error_callback=error_callback,
            workers=workers,
        )

        progress_bar.progress(1.0)
        elapsed = _time.time() - _load_start

        if not annotated_queries:
            status.update(label="No spectra could be loaded or annotated", state="error")
            return

        status_text.caption(f"Building DataFrames from {len(annotated_queries)} spectra...")

        peaks_lf = build_peaks_dataframe([(pid, df) for pid, _, df in annotated_queries])
        seq_lf = build_sequence_dataframe([(pid, s) for pid, s, _ in annotated_queries])
        heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated_queries])

        # Store everything in session state
        st.session_state[get_state_key('psm_lf')] = psm_lf
        st.session_state[get_state_key('psm_lookup')] = psm_lookup
        st.session_state[get_state_key('peaks_lf')] = peaks_lf
        st.session_state[get_state_key('seq_lf')] = seq_lf
        st.session_state[get_state_key('heatmap_lf')] = heatmap_lf

        st.session_state[get_state_key('annotated_queries')] = {
            pid: (spec, df) for pid, spec, df in annotated_queries
        }
        st.session_state[get_state_key('reference_store')] = reference_store

        import pandas as pd
        if comparison_rows:
            comparison_df = pd.DataFrame(comparison_rows)
        else:
            comparison_df = pd.DataFrame()
        st.session_state[get_state_key('comparison_df')] = comparison_df

        st.session_state[get_state_key('data_ready')] = True
        st.session_state[get_state_key('cache_version')] = \
            st.session_state.get(get_state_key('cache_version'), 0) + 1

        n_comparisons = len(comparison_rows)
        status.update(
            label=f"Loaded {len(annotated_queries)} spectra, {n_comparisons} comparisons "
                  f"in {format_time(elapsed)}",
            state="complete",
        )

    st.rerun()


def _get_koina_prediction(peptidoform: str, charge: int, annotation_kwargs: dict):
    """Generate Koina AI-predicted spectrum as a fallback reference.

    Returns (LoadedSpectrum-like, annotated_df) tuple or None.
    """
    try:
        from src.koina_prediction import predict_koina_spectrum, KOINA_AVAILABLE
        if not KOINA_AVAILABLE:
            return None

        import numpy as np
        from src.spectrum_loader import LoadedSpectrum
        from src.peak_annotation import prepare_spectrum_for_plot

        model_id = st.session_state.get(get_state_key("koina_model"), "ms2pip_timsTOF2024")
        ce = st.session_state.get(get_state_key("koina_ce"), 35)
        instrument = st.session_state.get(get_state_key("koina_instrument"), "TIMSTOF")

        pred_df = predict_koina_spectrum(
            peptidoform_str=peptidoform,
            charge=charge,
            model_id=model_id,
            collision_energy=ce,
            instrument_type=instrument,
        )

        if pred_df is None or pred_df.empty:
            return None

        # Extract base sequence (strip modifications)
        try:
            from psm_utils import Peptidoform as _PF
            base_seq = _PF(peptidoform.split("/")[0]).sequence
        except Exception:
            base_seq = peptidoform.split("/")[0]

        ref_spec = LoadedSpectrum(
            spectrum_id="koina_prediction",
            mz=pred_df["mz"].values.astype(np.float64),
            intensity=pred_df["intensity"].values.astype(np.float64),
            precursor_mz=0.0,
            precursor_charge=charge,
            peptidoform=peptidoform,
            base_sequence=base_seq,
            rt=0.0,
            source_file="Koina",
            is_reference=True,
        )

        ref_df = prepare_spectrum_for_plot(
            ref_spec.mz, ref_spec.intensity,
            ref_spec.peptidoform, ref_spec.precursor_charge,
            **annotation_kwargs,
        )

        return (ref_spec, ref_df)

    except Exception as e:
        logger.warning(f"Koina prediction failed for {peptidoform}: {e}")
        return None


def display_search_results(results: dict):
    """Display search results table."""
    from src.peptide_search import summarize_search_results
    summary_lf = summarize_search_results(results)
    summary_df = summary_lf.collect().to_pandas()
    if not summary_df.empty:
        display_df = summary_df.rename(columns={
            'base_sequence': 'Peptide', 'modification': 'Modification',
            'psm_count': 'PSMs', 'idxml_files': 'idXML Files',
            'mzml_resolved': 'mzML Resolved', 'unresolved': 'Unresolved'
        })
        display_df = display_df.drop(columns=['display_label'], errors='ignore')
        st.dataframe(display_df, width="stretch", hide_index=True)


def render_search_tab():
    idxmls = st.session_state.get(get_state_key('selected_idxmls'), [])
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    if not peptide_text:
        peptide_text = st.session_state.get(f"{STATE_PREFIX}_peptide_input", '')
    peptides = parse_peptide_list(peptide_text)

    if not idxmls:
        st.warning("Please select idXML files first (Setup tab)")
        return
    if not peptides:
        st.warning("Please enter target peptides first (Setup tab)")
        return

    if st.button("Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(idxmls, peptides)

    variant_groups = st.session_state.get(get_state_key('variant_groups'))
    if variant_groups:
        st.subheader("Search Results")
        display_search_results(variant_groups)

        total_matches = sum(len(vg.matches) for vg in variant_groups.values())
        st.success(f"Found {total_matches} PSMs across {len(variant_groups)} peptide variants. "
                   f"Click **Search & Load Spectra** below to load and annotate them.")
        mzmls = st.session_state.get(get_state_key('selected_mzmls'), [])
        if mzmls:
            st.subheader("Load Spectra & Compute Metrics")
            ref_sub = st.session_state.get(get_state_key('reference_substring'), '')
            if ref_sub:
                st.caption(f"Will pair query spectra with references (substring: '{ref_sub}') "
                           f"and compute 7 comparison metrics per pair.")
            else:
                koina_on = st.session_state.get(get_state_key('koina_enabled'), False)
                if koina_on:
                    st.caption("No reference configured -- will use Koina AI prediction as fallback.")
                else:
                    st.caption("No reference configured -- will annotate individual spectra.")
            if st.button("Search & Load Spectra", key=f"{STATE_PREFIX}_load_btn", type="primary"):
                run_search_and_load()

            if st.session_state.get(get_state_key('data_ready')):
                comparison_df = st.session_state.get(get_state_key('comparison_df'))
                if comparison_df is not None and not comparison_df.empty:
                    st.success(f"Data loaded with {len(comparison_df)} query/reference comparisons. "
                               f"Go to **Explore** tab.")
                else:
                    st.success("Data loaded. Go to **Explore** tab.")
        else:
            st.warning("Please select mzML files first (Setup tab)")
    else:
        st.caption("Search results will appear here after running the search.")


# ---------------------------------------------------------------------------
# Tab 3: Explore -- openms_insight + Plotly mirror plots + metrics
# ---------------------------------------------------------------------------

@st.fragment
def render_explore_tab():
    """Explore tab: interactive openms_insight components + mirror plots."""
    if not st.session_state.get(get_state_key('data_ready')):
        st.info("Run 'Search & Load' first to populate the explorer.")
        return

    explore_mode = st.radio(
        "View mode", options=["Interactive Explorer", "Mirror Plots", "Comparison Metrics"],
        horizontal=True, key=f"{STATE_PREFIX}_explore_mode"
    )

    if explore_mode == "Interactive Explorer":
        _render_openms_insight_explorer()
    elif explore_mode == "Mirror Plots":
        _render_mirror_plots()
    elif explore_mode == "Comparison Metrics":
        _render_comparison_metrics()


def _render_openms_insight_explorer():
    """openms_insight Table + LinePlot + SequenceView + Heatmap."""
    from openms_insight import Table, LinePlot, SequenceView, Heatmap, StateManager

    CACHE_PATH = st.session_state.get(get_state_key('export_dir'), '/tmp/msinsight_cache')
    os.makedirs(CACHE_PATH, exist_ok=True)

    cache_ver = st.session_state.get(get_state_key('cache_version'), 0)
    regen_key = get_state_key('cache_rendered_version')
    needs_regen = st.session_state.get(regen_key, -1) != cache_ver

    if 'msinsight_sm' not in st.session_state:
        st.session_state['msinsight_sm'] = StateManager(session_key="msinsight_state")
    sm = st.session_state['msinsight_sm']

    psm_lf = st.session_state.get(get_state_key('psm_lf'))
    peaks_lf = st.session_state.get(get_state_key('peaks_lf'))
    seq_lf = st.session_state.get(get_state_key('seq_lf'))
    heatmap_lf = st.session_state.get(get_state_key('heatmap_lf'))

    # Show currently selected PSM info
    current_state = sm.get_state_for_vue()
    selected_id = current_state.get("selected_psm")
    if selected_id is not None and psm_lf is not None:
        try:
            row = psm_lf.filter(pl.col("id") == selected_id).collect()
            if row.height > 0:
                label = row["display_label"][0]
                source = row["source_file"][0]
                st.info(f"Viewing: **{label}** (PSM {selected_id}) from {source}")
        except Exception:
            pass

    # PSM Table
    st.subheader("PSM Table")
    st.caption("Click a row to update the spectrum and sequence views below.")
    psm_columns = [
        {"field": "display_label", "title": "Peptide"},
        {"field": "charge", "title": "Z"},
        {"field": "score", "title": "Score", "formatter": "money", "formatterParams": {"precision": 4}},
        {"field": "source_file", "title": "Source File"},
    ]
    if needs_regen:
        table = Table(
            cache_id="psm_table",
            data=psm_lf,
            interactivity={"selected_psm": "id"},
            column_definitions=psm_columns,
            cache_path=CACHE_PATH,
            regenerate_cache=True,
            index_field="id",
        )
    else:
        table = Table(cache_id="psm_table", cache_path=CACHE_PATH)
    table(key="psm_table", state_manager=sm)

    # SequenceView
    st.subheader("Sequence Coverage")
    if needs_regen:
        peaks_for_sv = peaks_lf.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
        sv = SequenceView(
            cache_id="seq_view",
            sequence_data=seq_lf,
            peaks_data=peaks_for_sv,
            filters={"selected_psm": "psm_id"},
            annotation_config={
                "ion_types": st.session_state.get(get_state_key('ion_types'), ['b', 'y']),
                "tolerance": st.session_state.get(get_state_key('mz_tolerance'), 0.02),
                "tolerance_ppm": False,
            },
            cache_path=CACHE_PATH,
            height=350,
        )
    else:
        sv = SequenceView(cache_id="seq_view", cache_path=CACHE_PATH)
    sv(key="seq_view", state_manager=sm)

    # LinePlot linked to SequenceView for per-ion-type coloring
    st.subheader("Fragment Spectrum")
    if needs_regen:
        lineplot = LinePlot(
            cache_id="spectrum_plot",
            data=peaks_lf,
            filters={"selected_psm": "psm_id"},
            filter_defaults={"selected_psm": 0},
            cache_path=CACHE_PATH,
            regenerate_cache=True,
            x_column="x", y_column="y",
            highlight_column="highlight",
            annotation_column="annotation",
            x_label="m/z", y_label="Intensity",
        )
    else:
        lineplot = LinePlot(cache_id="spectrum_plot", cache_path=CACHE_PATH)
    lineplot(key="spectrum_plot", state_manager=sm, sequence_view_key="seq_view")

    # Heatmap
    st.subheader("RT x m/z Heatmap")
    if needs_regen:
        heatmap_comp = Heatmap(
            cache_id="rt_heatmap",
            data=heatmap_lf,
            x_column="x", y_column="y", intensity_column="intensity",
            x_label="Retention Time (s)", y_label="m/z",
            log_scale=True,
            cache_path=CACHE_PATH,
            regenerate_cache=True,
        )
    else:
        heatmap_comp = Heatmap(cache_id="rt_heatmap", cache_path=CACHE_PATH)
    heatmap_comp(key="rt_heatmap", state_manager=sm)

    if needs_regen:
        st.session_state[regen_key] = cache_ver


def _render_mirror_plots():
    """Interactive Plotly mirror plots for each PSM vs its references."""
    from src.mirror_plot import (
        create_plotly_mirror_plot, create_plotly_single_plot, ION_COLORS,
    )
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text

    annotated_queries = st.session_state.get(get_state_key('annotated_queries'), {})
    reference_store = st.session_state.get(get_state_key('reference_store'), {})
    variant_groups = st.session_state.get(get_state_key('variant_groups'), {})

    if not annotated_queries:
        st.info("No annotated spectra available.")
        return

    # Build list of PSMs for selection
    psm_options = []
    psm_id_to_vk = {}
    for vk, vg in variant_groups.items():
        for match in vg.matches:
            from src.spectrum_loader import spectrum_lookup_key
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for pid, (spec, _) in annotated_queries.items():
                if spectrum_lookup_key(spec) == key:
                    label = f"[{pid}] {vg.display_label} -- {Path(spec.source_file).stem} (RT: {spec.rt:.1f})"
                    psm_options.append((pid, label))
                    psm_id_to_vk[pid] = vk
                    break

    if not psm_options:
        for pid, (spec, _) in sorted(annotated_queries.items()):
            label = f"[{pid}] {spec.peptidoform} -- {Path(spec.source_file).stem}"
            psm_options.append((pid, label))

    selected_label = st.selectbox("Select PSM", options=[lbl for _, lbl in psm_options],
                                  key=f"{STATE_PREFIX}_mirror_psm_select")

    selected_pid = None
    for pid, lbl in psm_options:
        if lbl == selected_label:
            selected_pid = pid
            break

    if selected_pid is None or selected_pid not in annotated_queries:
        return

    query_spec, query_df = annotated_queries[selected_pid]
    vk = psm_id_to_vk.get(selected_pid)
    refs = reference_store.get(vk, []) if vk else []

    normalize = st.session_state.get(get_state_key('intensity_mode'), 'Relative') == 'Relative'

    if refs:
        for ref_spec, ref_df in refs:
            metrics = compute_all_metrics(query_df, ref_df)
            metrics_text = format_metrics_text(metrics)

            query_label = f"Query: {Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
            ref_source = Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina"
            ref_label = f"Ref: {ref_source} (RT: {ref_spec.rt:.1f})"
            title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"

            fig = create_plotly_mirror_plot(
                upper_spectrum=query_df,
                lower_spectrum=ref_df,
                title=title,
                upper_label=query_label,
                lower_label=ref_label,
                metrics_text=metrics_text,
                normalize=normalize,
            )
            st.plotly_chart(fig, width="stretch")

            with st.expander("Comparison Metrics", expanded=False):
                import pandas as pd
                from src.spectrum_metrics import METRIC_SHORT_TO_DISPLAY
                metric_data = {METRIC_SHORT_TO_DISPLAY[k]: f"{v:.4f}" if not (isinstance(v, float) and v != v) else "N/A"
                               for k, v in metrics.items() if k in METRIC_SHORT_TO_DISPLAY}
                metric_data["Shared Ions"] = str(metrics.get("shared", 0))
                metric_data["Annotated (Query)"] = str(metrics.get("ann_q", 0))
                metric_data["Annotated (Ref)"] = str(metrics.get("ann_r", 0))
                st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)
    else:
        title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
        label = f"{Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
        fig = create_plotly_single_plot(
            spectrum=query_df, title=title, label=label, normalize=normalize,
        )
        st.plotly_chart(fig, width="stretch")


def _render_comparison_metrics():
    """Show comparison metrics table and summary statistics."""
    import pandas as pd

    comparison_df = st.session_state.get(get_state_key('comparison_df'))
    if comparison_df is None or comparison_df.empty:
        st.info("No comparison metrics available. Configure a reference substring and re-run.")
        return

    st.subheader("Comparison Summary")
    st.caption(f"{len(comparison_df)} query/reference comparisons")

    # Display columns
    display_cols = ["Peptide", "Modification", "Charge", "Query", "Reference",
                    "R\u00b2 match", "R\u00b2 binned", "Cosine", "Entropy",
                    "Dice", "Jaccard", "SharedRatio", "shared", "ann_q", "ann_r"]
    available_cols = [c for c in display_cols if c in comparison_df.columns]
    st.dataframe(comparison_df[available_cols], width="stretch", hide_index=True)

    # Ion match statistics (if available)
    ims_cols = [c for c in comparison_df.columns if c.startswith("ims_")]
    if ims_cols:
        with st.expander(f"Ion Match Statistics ({len(ims_cols)} metrics)", expanded=False):
            ims_display_cols = ["Peptide", "Query", "Reference"] + ims_cols
            ims_available = [c for c in ims_display_cols if c in comparison_df.columns]
            st.dataframe(comparison_df[ims_available], width="stretch", hide_index=True)

    # Summary statistics per peptidoform
    st.markdown("---")
    st.subheader("Summary Statistics")
    from src.spectrum_metrics import METRIC_SHORT_TO_DISPLAY
    metric_cols = [v for v in METRIC_SHORT_TO_DISPLAY.values() if v in comparison_df.columns]
    if metric_cols:
        summary = comparison_df.groupby(["Peptide", "Modification", "Charge"])[metric_cols].agg(
            ["count", "median", "mean", "std"]
        )
        summary.columns = [f"{m} {s}" for m, s in summary.columns]
        st.dataframe(summary, width="stretch")


# ---------------------------------------------------------------------------
# Tab 4: Export (with parallel generation, versioning, history, gallery)
# ---------------------------------------------------------------------------

def render_export_tab():
    """Export: Parquet, CSV/XLSX, mirror plots, summary plots, HTML, ZIP."""
    if not st.session_state.get(get_state_key('data_ready')):
        st.info("Load data first (Search & Load tab).")
        return

    _ws = st.session_state.get('workspace', None)
    _export_fallback = str(Path(_ws) / "exports") if _ws else "./exports"
    export_dir = st.text_input(
        "Export Directory",
        value=st.session_state.get(get_state_key('export_dir'), _export_fallback),
        key=f"{STATE_PREFIX}_export_dir_input"
    )
    st.session_state[get_state_key('export_dir')] = export_dir

    export_col1, export_col2 = st.columns(2)
    with export_col1:
        export_parquet = st.checkbox("Parquet DataFrames", value=True,
                                     key=f"{STATE_PREFIX}_exp_parquet")
        export_csv = st.checkbox("Comparison Summary (CSV/XLSX)", value=True,
                                 key=f"{STATE_PREFIX}_exp_csv")
        export_html = st.checkbox("Interactive Plots (HTML)", value=False,
                                  key=f"{STATE_PREFIX}_exp_html")
    with export_col2:
        export_mirror = st.checkbox("Mirror Plots (PNG/SVG)", value=True,
                                    key=f"{STATE_PREFIX}_exp_mirror")
        export_summary = st.checkbox("Summary Analysis Plots", value=True,
                                     key=f"{STATE_PREFIX}_exp_summary")

    if st.button("Export All", key=f"{STATE_PREFIX}_export_btn", type="primary"):
        _run_export(export_dir, export_parquet, export_csv, export_mirror,
                    export_summary, export_html)

    # Show export history
    st.markdown("---")
    _render_export_history(export_dir)

    # ZIP download (if export dir exists and has files)
    export_path = Path(export_dir)
    if export_path.exists() and any(export_path.rglob("*")):
        st.markdown("---")
        st.subheader("Download")
        _render_zip_download(export_path)


def _build_run_settings() -> dict:
    """Capture current settings for run_settings.json."""
    peptide_text = st.session_state.get(get_state_key('target_peptides'), '')
    peptides = parse_peptide_list(peptide_text)
    return {
        "idxml_files": st.session_state.get(get_state_key('selected_idxmls'), []),
        "mzml_files": st.session_state.get(get_state_key('selected_mzmls'), []),
        "peptides": peptides,
        "reference_substring": st.session_state.get(get_state_key('reference_substring'), ''),
        "ion_types": st.session_state.get(get_state_key('ion_types'), ['b', 'y']),
        "mz_tolerance": st.session_state.get(get_state_key('mz_tolerance'), 0.02),
        "add_losses": st.session_state.get(get_state_key('add_losses'), False),
        "add_first_prefix_ion": st.session_state.get(get_state_key('add_first_prefix_ion'), True),
        "remove_precursor_peaks": st.session_state.get(get_state_key('remove_precursor_peaks'), False),
        "intensity_mode": st.session_state.get(get_state_key('intensity_mode'), 'Relative'),
        "koina_enabled": st.session_state.get(get_state_key('koina_enabled'), False),
        "_created_at": datetime.now().isoformat(),
        "_last_run": datetime.now().isoformat(),
    }


def _run_export(export_dir: str, do_parquet: bool, do_csv: bool,
                do_mirror: bool, do_summary: bool, do_html: bool):
    """Execute all selected exports with versioned directory and parallel generation."""
    import pandas as pd
    import traceback as _tb

    # Create versioned export directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    export_root = Path(export_dir) / f"export_{timestamp}"
    export_root.mkdir(parents=True, exist_ok=True)

    saved_files = []

    with st.status("Exporting...", expanded=True) as status:
        # Save run_settings.json
        run_settings = _build_run_settings()
        settings_path = export_root / "run_settings.json"
        with open(settings_path, "w") as f:
            json.dump(run_settings, f, indent=2, default=str)
        saved_files.append(str(settings_path))

        # 1. Parquet DataFrames
        if do_parquet:
            st.caption("Saving Parquet DataFrames...")
            for name in ('psm_lf', 'peaks_lf', 'seq_lf', 'heatmap_lf'):
                lf = st.session_state.get(get_state_key(name))
                if lf is not None:
                    out_file = export_root / f"{name.replace('_lf', '')}.parquet"
                    lf.collect().write_parquet(out_file)
                    saved_files.append(str(out_file))

        # 2. Comparison summary CSV/XLSX
        if do_csv:
            comparison_df = st.session_state.get(get_state_key('comparison_df'))
            if comparison_df is not None and not comparison_df.empty:
                st.caption("Saving comparison summary...")
                csv_path = export_root / "comparison_summary.csv"
                comparison_df.to_csv(csv_path, index=False)
                saved_files.append(str(csv_path))

                try:
                    xlsx_path = export_root / "comparison_summary.xlsx"
                    comparison_df.to_excel(xlsx_path, index=False, engine="openpyxl")
                    saved_files.append(str(xlsx_path))
                except ImportError:
                    logger.warning("openpyxl not installed -- skipping XLSX export")
                except Exception as e:
                    logger.warning(f"XLSX export failed: {e}")

        # 3. Mirror plots
        if do_mirror or do_html:
            st.caption("Generating mirror plots...")
            try:
                n_plots, plot_files = _export_mirror_plots(export_root, do_html)
                saved_files.extend(plot_files)
                st.caption(f"  Generated {n_plots} mirror plot file(s)")
            except Exception as e:
                st.error(f"Mirror plot export failed: {e}")
                logger.error(f"Mirror plot export failed:\n{_tb.format_exc()}")
        else:
            st.caption("Mirror plots: skipped (unchecked)")

        # 4. Summary analysis plots
        if do_summary:
            comparison_df = st.session_state.get(get_state_key('comparison_df'))
            if comparison_df is not None and not comparison_df.empty:
                st.caption("Generating summary plots...")
                from src.summary_plots import generate_all_summary_plots
                try:
                    results = generate_all_summary_plots(
                        source=comparison_df,
                        output_dir=str(export_root),
                        formats=["png", "svg"],
                    )
                    n_summary = sum(len(v) for v in results.values())
                    st.caption(f"  Generated {n_summary} summary plot(s)")
                    for pt, paths in results.items():
                        for fmt, p in paths.items():
                            saved_files.append(p)
                except Exception as e:
                    logger.error(f"Summary plot generation failed: {e}")
                    st.warning(f"Summary plots failed: {e}")

        status.update(label=f"Exported {len(saved_files)} files to {export_root}", state="complete")

    if saved_files:
        with st.expander(f"Exported files ({len(saved_files)})", expanded=False):
            for f in saved_files:
                st.caption(f"  {f}")

    # Show plot gallery for this export
    _render_plot_gallery(export_root)


def _export_mirror_plots(export_root: Path, do_html: bool):
    """Generate mirror/single plots organized by peptidoform.

    Returns (file_count, list_of_file_paths).  All errors are surfaced
    via ``st.warning`` so nothing fails silently.
    """
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text
    from src.summary_plots import sanitize_filename
    from src.mirror_plot import (
        create_mirror_plot, create_single_spectrum_plot,
        figure_to_bytes, create_plotly_mirror_plot,
    )
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    annotated_queries = st.session_state.get(get_state_key('annotated_queries'), {})
    reference_store = st.session_state.get(get_state_key('reference_store'), {})
    variant_groups = st.session_state.get(get_state_key('variant_groups'), {})

    if not annotated_queries:
        st.warning("No annotated spectra in session -- re-run **Search & Load** first.")
        return 0, []

    st.caption(f"  {len(annotated_queries)} annotated spectra in session")

    normalize = st.session_state.get(get_state_key('intensity_mode'), 'Relative') == 'Relative'

    png_dir = export_root / "png"
    svg_dir = export_root / "svg"
    html_dir = export_root / "html" if do_html else None

    # Build psm_id -> variant_key mapping
    from src.spectrum_loader import spectrum_lookup_key
    psm_id_to_vk = {}
    for vk, vg in variant_groups.items():
        for match in vg.matches:
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for pid, (spec, _) in annotated_queries.items():
                if spectrum_lookup_key(spec) == key:
                    psm_id_to_vk[pid] = vk
                    break

    saved_files = []
    errors = []

    for pid, (query_spec, query_df) in sorted(annotated_queries.items()):
        vk = psm_id_to_vk.get(pid)
        refs = reference_store.get(vk, []) if vk else []
        folder_name = sanitize_filename(query_spec.peptidoform)

        if refs:
            for ref_spec, ref_df in refs:
                try:
                    metrics = compute_all_metrics(query_df, ref_df)
                    metrics_text = format_metrics_text(metrics)
                    ref_source = Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina"
                    safe_name = sanitize_filename(
                        f"{Path(query_spec.source_file).stem}_vs_{ref_source}"
                    )
                    title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
                    query_label = f"Query: {Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
                    ref_label = f"Ref: {ref_source} (RT: {ref_spec.rt:.1f})"

                    fig = create_mirror_plot(
                        upper_spectrum=query_df, lower_spectrum=ref_df,
                        title=title, upper_label=query_label, lower_label=ref_label,
                        metrics_text=metrics_text, normalize=normalize,
                    )

                    for fmt_dir, fmt in [(png_dir, "png"), (svg_dir, "svg")]:
                        out_dir = fmt_dir / folder_name
                        out_dir.mkdir(parents=True, exist_ok=True)
                        out_file = out_dir / f"{safe_name}.{fmt}"
                        out_file.write_bytes(figure_to_bytes(fig, format=fmt))
                        saved_files.append(str(out_file))

                    plt.close(fig)

                    if do_html and html_dir:
                        plotly_fig = create_plotly_mirror_plot(
                            upper_spectrum=query_df, lower_spectrum=ref_df,
                            title=title, upper_label=query_label,
                            lower_label=ref_label, metrics_text=metrics_text,
                            normalize=normalize,
                        )
                        out_dir = html_dir / folder_name
                        out_dir.mkdir(parents=True, exist_ok=True)
                        html_path = out_dir / f"{safe_name}.html"
                        plotly_fig.write_html(str(html_path), include_plotlyjs="cdn")
                        saved_files.append(str(html_path))

                except Exception as e:
                    errors.append(f"{folder_name}/{safe_name}: {e}")
                    logger.error(f"Mirror plot failed for PSM {pid}: {e}")
        else:
            try:
                safe_name = sanitize_filename(Path(query_spec.source_file).stem)
                title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
                label = f"{Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"

                fig = create_single_spectrum_plot(
                    spectrum=query_df, title=title, label=label, normalize=normalize,
                )

                for fmt_dir, fmt in [(png_dir, "png"), (svg_dir, "svg")]:
                    out_dir = fmt_dir / folder_name
                    out_dir.mkdir(parents=True, exist_ok=True)
                    out_file = out_dir / f"{safe_name}.{fmt}"
                    out_file.write_bytes(figure_to_bytes(fig, format=fmt))
                    saved_files.append(str(out_file))

                plt.close(fig)

            except Exception as e:
                errors.append(f"{folder_name}/{safe_name}: {e}")
                logger.error(f"Single plot failed for PSM {pid}: {e}")

    if errors:
        st.warning(f"{len(errors)} plot(s) failed:")
        for err in errors[:10]:
            st.caption(f"  {err}")

    return len(saved_files), saved_files


def _render_plot_gallery(root_dir: Path) -> None:
    """Display plots saved in root_dir organized by peptidoform then reference."""
    format_dirs = [(root_dir / "svg", "svg"), (root_dir / "png", "png")]
    existing_format_dirs = [(d, fmt) for d, fmt in format_dirs if d.exists()]

    if not existing_format_dirs:
        return

    browse_dir, browse_fmt = existing_format_dirs[0]
    html_base_dir = root_dir / "html"

    peptide_dirs = sorted([d for d in browse_dir.iterdir() if d.is_dir()])
    if not peptide_dirs:
        return

    st.markdown("---")
    st.subheader("Browse Plots")

    total_plots = sum(len(list(d.glob(f"*.{browse_fmt}"))) for d in peptide_dirs)
    st.caption(f"{len(peptide_dirs)} peptide variants, {total_plots} plots")

    for peptide_dir in peptide_dirs:
        plot_files = sorted(peptide_dir.glob(f"*.{browse_fmt}"))
        if not plot_files:
            continue

        with st.expander(f"{peptide_dir.name} ({len(plot_files)} plots)", expanded=False):
            ref_groups: dict = {}
            for pf in plot_files:
                ref_stem = pf.stem.rsplit("_vs_", 1)[1] if "_vs_" in pf.stem else "(no reference)"
                ref_groups.setdefault(ref_stem, []).append(pf)

            for ref_stem, ref_files in sorted(ref_groups.items()):
                sub_label = (
                    f"Reference: {ref_stem}" if ref_stem != "(no reference)"
                    else "Single spectra"
                )
                with st.expander(f"{sub_label} ({len(ref_files)})", expanded=False):
                    for plot_file in ref_files[:20]:
                        # Show the image
                        if browse_fmt == "png":
                            st.image(str(plot_file), caption=plot_file.stem, width="stretch")
                        elif browse_fmt == "svg":
                            try:
                                svg_content = plot_file.read_text()
                                st.markdown(svg_content, unsafe_allow_html=True)
                            except Exception:
                                st.caption(f"Could not render: {plot_file.name}")

                        # Download buttons
                        run_id = root_dir.name
                        dl_col1, dl_col2, dl_col3 = st.columns(3)
                        with dl_col1:
                            png_path = root_dir / "png" / peptide_dir.name / f"{plot_file.stem}.png"
                            if png_path.exists():
                                st.download_button(
                                    "PNG", data=png_path.read_bytes(),
                                    file_name=f"{plot_file.stem}.png",
                                    mime="image/png",
                                    key=f"dl_png_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )
                        with dl_col2:
                            svg_path = root_dir / "svg" / peptide_dir.name / f"{plot_file.stem}.svg"
                            if svg_path.exists():
                                st.download_button(
                                    "SVG", data=svg_path.read_bytes(),
                                    file_name=f"{plot_file.stem}.svg",
                                    mime="image/svg+xml",
                                    key=f"dl_svg_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )
                        with dl_col3:
                            html_file = html_base_dir / peptide_dir.name / f"{plot_file.stem}.html"
                            if html_file.exists():
                                st.download_button(
                                    "HTML", data=html_file.read_bytes(),
                                    file_name=f"{plot_file.stem}.html",
                                    mime="text/html",
                                    key=f"dl_html_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )

                    if len(ref_files) > 20:
                        st.caption(f"... and {len(ref_files) - 20} more")


def _render_export_history(export_dir: str):
    """Browse previous export runs, reload their settings, or delete them."""
    st.subheader("Export History")

    export_path = Path(export_dir)
    if not export_path.exists():
        st.caption(f"No exports directory found at `{export_dir}`.")
        return

    runs = []
    for settings_file in sorted(export_path.glob("*/run_settings.json"), reverse=True):
        try:
            with open(settings_file) as f:
                run_cfg = json.load(f)
            runs.append((settings_file.parent, run_cfg))
        except Exception:
            continue

    if not runs:
        st.caption(f"No completed export runs found in `{export_dir}`.")
        return

    st.caption(f"{len(runs)} run(s) in `{export_dir}`")

    for run_dir, run_cfg in runs:
        peptides = run_cfg.get('peptides', [])
        idxml_count = len(run_cfg.get('idxml_files', []))
        mzml_count = len(run_cfg.get('mzml_files', []))
        created = (run_cfg.get('_created_at') or '')[:10]

        pep_preview = ', '.join(peptides[:3])
        if len(peptides) > 3:
            pep_preview += f' ... +{len(peptides) - 3}'

        label = (
            f"**{run_dir.name}**"
            f" -- {created}"
            f" -- {len(peptides)} peptide(s)"
            f" -- {idxml_count} idXML / {mzml_count} mzML"
        )

        with st.expander(label, expanded=False):
            col_meta, col_actions = st.columns([3, 1])

            with col_meta:
                st.markdown(f"**Peptides:** {pep_preview or '--'}")
                st.markdown(f"**Files:** {idxml_count} idXML, {mzml_count} mzML")
                ref_sub = run_cfg.get('reference_substring') or '--'
                st.markdown(f"**Reference substring:** `{ref_sub}`")
                st.markdown(f"**Ion types:** {', '.join(run_cfg.get('ion_types', []))}")
                st.markdown(f"**m/z tolerance:** {run_cfg.get('mz_tolerance', '--')} Da")
                created_full = (run_cfg.get('_created_at') or '--')[:19]
                last_run_full = (run_cfg.get('_last_run') or '--')[:19]
                st.caption(f"Created: {created_full}  |  Last run: {last_run_full}")

            with col_actions:
                if st.button("Reload settings", key=f"reload_{run_dir.name}"):
                    _apply_run_settings(run_cfg)
                    st.success("Settings loaded -- switch to Setup tab to review.")
                    st.rerun()

                confirm_key = f"confirm_delete_{run_dir.name}"
                if st.session_state.get(confirm_key):
                    st.warning("Delete all files in this run?")
                    c1, c2 = st.columns(2)
                    with c1:
                        if st.button("Yes, delete", key=f"yes_del_{run_dir.name}", type="primary"):
                            shutil.rmtree(run_dir)
                            del st.session_state[confirm_key]
                            st.rerun()
                    with c2:
                        if st.button("Cancel", key=f"cancel_del_{run_dir.name}"):
                            del st.session_state[confirm_key]
                            st.rerun()
                else:
                    if st.button("Delete run", key=f"delete_{run_dir.name}"):
                        st.session_state[confirm_key] = True
                        st.rerun()

            # Plot gallery for this run
            _render_plot_gallery(run_dir)


def _apply_run_settings(run_cfg: dict):
    """Load a run_settings.json dict back into session state."""
    mapping = {
        'selected_idxmls':        run_cfg.get('idxml_files', []),
        'selected_mzmls':         run_cfg.get('mzml_files', []),
        'target_peptides':        '\n'.join(run_cfg.get('peptides', [])),
        'reference_substring':    run_cfg.get('reference_substring', ''),
        'ion_types':              run_cfg.get('ion_types', ['b', 'y']),
        'mz_tolerance':           run_cfg.get('mz_tolerance', 0.02),
        'add_losses':             run_cfg.get('add_losses', False),
        'add_first_prefix_ion':   run_cfg.get('add_first_prefix_ion', True),
        'remove_precursor_peaks': run_cfg.get('remove_precursor_peaks', False),
        'intensity_mode':         run_cfg.get('intensity_mode', 'Relative'),
        'koina_enabled':          run_cfg.get('koina_enabled', False),
    }
    for key, value in mapping.items():
        st.session_state[get_state_key(key)] = value


@st.fragment
def _render_zip_download(export_root: Path) -> None:
    """Deferred zip download of the export directory."""
    def _make_zip() -> bytes:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for file in export_root.rglob("*"):
                if not file.is_file():
                    continue
                zf.write(file, file.relative_to(export_root.parent))
        return buf.getvalue()

    st.download_button(
        label="Download Export as ZIP",
        data=_make_zip,
        file_name=f"{export_root.name}.zip",
        mime="application/zip",
    )


# ---------------------------------------------------------------------------
# Main page layout
# ---------------------------------------------------------------------------

def render():
    """Main render function for the page."""
    init_session_state(params)

    st.title("Peptide Mirror Analysis")
    st.caption("Compare experimental spectra against reference spectra with interactive visualizations")

    tabs = st.tabs(["Setup", "Search & Load", "Explore", "Export"])

    with tabs[0]:
        render_setup_tab()
    with tabs[1]:
        render_search_tab()
    with tabs[2]:
        render_explore_tab()
    with tabs[3]:
        render_export_tab()


# Run render when page is loaded
render()
