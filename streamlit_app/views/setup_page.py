"""Setup tab — pick idXML/mzML files, enter peptides, configure reference + annotation + Koina.

Registered in app.py via st.Page(); outputs land in params and feed Load & Annotate."""
from pathlib import Path

import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from streamlit_app.views.widgets.file_browser import render_directory_browser
from streamlit_app.views.widgets.peptide_input import render_peptide_input
from core.defaults import DEFAULTS, INSTRUMENT_TYPES, KOINA_MODELS
from core.helpers import (
    DuplicateBasenameError,
    available_cpu_count,
    basenames_shared_across_dirs,
    ensure_unique_basenames,
)
from core.spectra import get_koina_settings


def _get_nested_value(source: dict, *keys, default=None):
    """Walk nested dicts safely — _get_nested_value(d, 'a', 'b') returns d['a']['b'] or *default*."""
    current = source
    for key in keys:
        if isinstance(current, dict):
            current = current.get(key)
        else:
            return default
    return current if current is not None else default


params = page_setup("Setup")


def render_file_selection():
    """Render side-by-side directory browsers for idXML (identifications) and mzML (spectra)."""
    # Resolve starting directories: workspace override beats settings.json default
    _settings = st.session_state.get("settings", {})
    data_dir_default = _get_nested_value(_settings, "data_directories", "default", default="/")

    idxml_dir = params.get('idxml_browser_dir', data_dir_default)
    mzml_dir = params.get('mzml_browser_dir', data_dir_default)

    # Render: two columns side-by-side, idXML left and mzML right
    st.header("File Selection")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_idxml",
                                 extensions=['.idXML', '.idxml'],
                                 state_key='selected_idxmls',
                                 default_dir=idxml_dir, params=params)
    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_mzml",
                                 extensions=['.mzML', '.mzml'],
                                 state_key='selected_mzmls',
                                 default_dir=mzml_dir, params=params)

    # Persist browser directories back into params and session state for next rerun
    for dir_param, dir_ss_key, fallback in [
        ('idxml_browser_dir', f"{STATE_PREFIX}_idxml_dir", idxml_dir),
        ('mzml_browser_dir', f"{STATE_PREFIX}_mzml_dir", mzml_dir),
    ]:
        selected_directory = st.session_state.get(dir_ss_key, fallback)
        params[dir_param] = selected_directory
        st.session_state[dir_param] = selected_directory

    # Validate: refuse only the exact same file selected twice. Same-basename files
    # from different source folders are allowed — downstream resolution and provenance
    # disambiguate them by source directory — but surface an info note so the user
    # knows the folder is what distinguishes them.
    for kind, state_key in [("idXML", 'selected_idxmls'), ("mzML", 'selected_mzmls')]:
        selected = st.session_state.get(state_key, [])
        try:
            ensure_unique_basenames(selected, kind)
        except DuplicateBasenameError as error:
            st.error(str(error))
        shared = basenames_shared_across_dirs(selected)
        if shared:
            st.info(
                f"{len(shared)} {kind} basename(s) appear in more than one folder "
                f"(e.g. **{next(iter(shared))}**). These will be kept distinct by their "
                f"source folder, so identifications from each can be compared."
            )

    save_params(params)


def render_reference_selection():
    """Configure the substring that splits selected mzMLs into reference vs query files."""
    # Guard: require mzML selections before showing the reference UI
    st.subheader("Reference Identification")
    mzmls = st.session_state.get('selected_mzmls', params.get('selected_mzmls', []))
    if not mzmls:
        st.warning("Please select mzML files first")
        return

    # Substring input: persist across page navigation via a non-widget key
    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")
    persisted_key = get_state_key("reference_substring_persisted")
    initial_value = st.session_state.get(persisted_key, params.get('reference_substring', ''))
    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=initial_value,
        key="reference_substring", placeholder="e.g., Synpep"
    )
    # Mirror into a non-widget key — Streamlit garbage-collects widget keys
    # when the owning page unmounts, which would lose the value on navigation.
    st.session_state[persisted_key] = ref_substring
    # Seed into params so save_params persists even in legacy workspaces whose
    # params.json predates the snake_case key migration.
    params['reference_substring'] = ref_substring

    # Preview: show the reference/query split (or an info note when blank).
    # Match against the full path — mirroring the loader's _split_reference_and_query_keys
    # — so a substring like "illumina" that lives in the folder (not the basename) still
    # splits correctly; labels show <folder>/<file> since basenames may repeat across folders.
    if ref_substring:
        needle = ref_substring.lower()
        labelled = [(f"{Path(m).parent.name}/{Path(m).name}", needle in str(m).lower()) for m in mzmls]
        reference_files = [label for label, is_ref in labelled if is_ref]
        query_files = [label for label, is_ref in labelled if not is_ref]

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
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")
    else:
        st.info("No reference identifier -- individual annotated spectrum plots only.")

    save_params(params)


def render_settings():
    """Render annotation, display, performance, and Koina prediction widgets — all writing into params."""
    st.subheader("Settings")
    available_ions = _get_nested_value(st.session_state.get("settings", {}), "ion_types", "available",
                                 default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

    # Ion types + m/z tolerance: two-column row at the top
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Ion Types**")
        current_ions = params.get('ion_types', ['b', 'y'])
        st.multiselect("Select ion types", options=available_ions,
                       default=[i for i in current_ions if i in available_ions],
                       key="ion_types")
    with col2:
        st.markdown("**m/z Tolerance**")
        st.number_input("Tolerance (Da)", min_value=0.001, max_value=1.0,
                        value=float(params.get('mz_tolerance', DEFAULTS.mz_tolerance)),
                        step=0.005, format="%.3f", key="mz_tolerance")

    # Annotation: neutral-loss toggle + raw-vs-relative intensity mode
    st.markdown("**Annotation Settings**")
    st.checkbox("Add neutral losses (H2O, NH3)",
                value=params.get('add_neutral_losses', DEFAULTS.add_neutral_losses),
                key="add_neutral_losses")
    st.radio("**Intensity Display Mode**",
             options=["Relative", "Raw"],
             captions=["Relative to the highest intensity", "Raw intensity"],
             index=0 if params.get('intensity_mode', DEFAULTS.intensity_mode) == 'Relative' else 1,
             horizontal=True, key="intensity_mode", label_visibility="visible")

    # Performance: ProcessPoolExecutor worker count for parallel pipeline stages.
    # Cap the slider at the cores this host can actually use, and clamp any saved
    # value into range so a workspace carried over from a bigger machine still loads.
    st.markdown("**Performance**")
    max_workers = available_cpu_count()
    if max_workers <= 1:
        # Single-core host: a slider needs min < max, so pin workers and say why.
        st.session_state["workers"] = 1
        st.caption("Parallel workers: 1 (only one core available on this host)")
    else:
        workers_value = min(params.get('workers', DEFAULTS.workers), max_workers)
        st.slider("Parallel workers", min_value=1, max_value=max_workers,
                  value=workers_value, key="workers",
                  help=f"Number of parallel processes for plot generation "
                       f"(~30 MB per worker; {max_workers} cores available)")

    # Koina AI: enable toggle reveals model / collision-energy / instrument widgets
    st.markdown("---")
    st.markdown("**Koina AI Prediction**")
    st.caption(
        "When enabled, spectra without an experimental reference are compared "
        "against AI-predicted fragmentation patterns (intensity + m/z)."
    )
    koina_defaults = get_koina_settings(st.session_state.get("settings", {}))
    koina_enabled = st.checkbox(
        "Enable Koina fallback reference",
        value=params.get("koina_enabled", koina_defaults["enabled"]),
        key="koina_enabled",
        help="Use AI-predicted spectrum as reference when no experimental reference mzML exists",
    )
    if koina_enabled:
        # Model / collision energy / instrument: three-column row of selectors
        k_col1, k_col2, k_col3 = st.columns(3)
        with k_col1:
            model_names = list(KOINA_MODELS.keys())
            model_ids = list(KOINA_MODELS.values())
            current_model_id = params.get("koina_model", koina_defaults["model"])
            default_idx = model_ids.index(current_model_id) if current_model_id in model_ids else 0
            selected_model_name = st.selectbox(
                "Prediction Model", options=model_names, index=default_idx,
                key=f"{STATE_PREFIX}_koina_model_select",
                help="Deep-learning model for fragment ion prediction",
            )
            st.session_state["koina_model"] = KOINA_MODELS[selected_model_name]
        with k_col2:
            st.number_input(
                "Collision Energy (eV)", min_value=1, max_value=100,
                value=params.get("koina_collision_energy", koina_defaults["collision_energy"]),
                key="koina_collision_energy",
            )
        with k_col3:
            current_inst = params.get("koina_instrument_type", koina_defaults["instrument_type"])
            inst_idx = INSTRUMENT_TYPES.index(current_inst) if current_inst in INSTRUMENT_TYPES else 0
            st.selectbox(
                "Instrument Type", options=INSTRUMENT_TYPES, index=inst_idx,
                key="koina_instrument_type",
            )

    save_params(params)


def render():
    """Render the four Setup sections top-to-bottom — files, peptides, reference, settings."""
    # Restore analysis: rehydrate workspace state if Load & Annotate was already run
    restore_analysis_if_needed()
    st.title("Setup")
    st.caption("Configure files, peptides, reference, and annotation settings")

    # Files: idXML + mzML directory browsers
    render_file_selection()
    st.markdown("---")
    # Peptides: target sequence list input (handled by shared widget)
    render_peptide_input(params)
    st.markdown("---")
    # Reference: substring that splits mzMLs into reference vs query
    render_reference_selection()
    st.markdown("---")
    # Settings: annotation, performance, Koina
    render_settings()

    save_params(params)


render()
