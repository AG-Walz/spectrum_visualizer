"""
Setup Page — File selection, peptide input, reference config, annotation settings.
"""
from pathlib import Path

import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from streamlit_app.views.widgets.file_browser import render_directory_browser
from streamlit_app.views.widgets.peptide_input import render_peptide_input
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.text_utils import (
    DuplicateBasenameError,
    ensure_unique_basenames,
    get_nested,
)

try:
    from mspecviz.spectra import (
        KOINA_AVAILABLE, KOINA_MODELS, INSTRUMENT_TYPES,
        get_koina_settings,
    )
except ImportError:
    KOINA_AVAILABLE = False
    KOINA_MODELS = {}
    INSTRUMENT_TYPES = []

    def get_koina_settings(s):
        return {"enabled": False, "model": DEFAULTS.koina_model,
                "collision_energy": DEFAULTS.koina_collision_energy,
                "instrument_type": DEFAULTS.koina_instrument_type}

params = page_setup()


def file_selection_fragment():
    """File browser for idXML and mzML files."""
    _settings = st.session_state.get("settings", {})
    data_dir_default = get_nested(_settings, "data_directories", "default", default="/")

    idxml_dir = params.get('idxml_browser_dir', data_dir_default)
    mzml_dir = params.get('mzml_browser_dir', data_dir_default)

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

    for dir_param, dir_ss_key, fallback in [
        ('idxml_browser_dir', f"{STATE_PREFIX}_idxml_dir", idxml_dir),
        ('mzml_browser_dir', f"{STATE_PREFIX}_mzml_dir", mzml_dir),
    ]:
        selected_directory = st.session_state.get(dir_ss_key, fallback)
        params[dir_param] = selected_directory
        st.session_state[dir_param] = selected_directory

    for kind, state_key in [("idXML", 'selected_idxmls'), ("mzML", 'selected_mzmls')]:
        try:
            ensure_unique_basenames(st.session_state.get(state_key, []), kind)
        except DuplicateBasenameError as exc:
            st.error(str(exc))

    save_params(params)


def reference_selection_fragment():
    """Reference substring configuration."""
    st.subheader("Reference Identification")
    mzmls = st.session_state.get('selected_mzmls', params.get('selected_mzmls', []))
    if not mzmls:
        st.warning("Please select mzML files first")
        return

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")
    persisted_key = get_state_key("reference_substring_persisted")
    initial_value = st.session_state.get(persisted_key, params.get('reference_substring', ''))
    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=initial_value,
        key="reference_substring", placeholder="e.g., Synpep"
    )
    # Mirror into a non-widget key so the value survives page navigation
    # (Streamlit garbage-collects widget keys when the owning page unmounts).
    st.session_state[persisted_key] = ref_substring
    # Seed into params so save_params picks it up even in legacy workspaces
    # whose params.json predates the snake_case key migration.
    params['reference_substring'] = ref_substring

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
        elif not reference_files:
            st.warning("No files match the reference identifier")
        elif not query_files:
            st.warning("All files match as reference - no query files")
    else:
        st.info("No reference identifier -- individual annotated spectrum plots only.")

    save_params(params)


def settings_fragment():
    """Annotation, display, performance, and Koina settings fragment."""
    st.subheader("Settings")
    available_ions = get_nested(st.session_state.get("settings", {}), "ion_types", "available",
                                default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

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

    st.markdown("**Annotation Settings**")
    st.checkbox("Add neutral losses (H2O, NH3)",
                value=params.get('add_neutral_losses', DEFAULTS.add_neutral_losses),
                key="add_neutral_losses")

    st.radio("**Intensity Display Mode**",
             options=["Relative", "Raw"],
             captions=["Relative to the highest intensity", "Raw intensity"],
             index=0 if params.get('intensity_mode', DEFAULTS.intensity_mode) == 'Relative' else 1,
             horizontal=True, key="intensity_mode", label_visibility="visible")

    st.markdown("**Performance**")
    st.slider("Parallel workers", min_value=1, max_value=32,
              value=params.get('workers', DEFAULTS.workers), key="workers",
              help="Number of parallel processes for plot generation (~30 MB per worker)")

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
        disabled=not KOINA_AVAILABLE,
    )

    if not KOINA_AVAILABLE:
        st.warning("koinapy is not installed. Install it to enable Koina predictions.")

    if koina_enabled and KOINA_AVAILABLE:
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
    restore_analysis_if_needed()
    st.title("Setup")
    st.caption("Configure files, peptides, reference, and annotation settings")

    file_selection_fragment()
    st.markdown("---")
    render_peptide_input(params)
    st.markdown("---")
    reference_selection_fragment()
    st.markdown("---")
    settings_fragment()

    save_params(params)


render()
