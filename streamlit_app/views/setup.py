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
from mspecviz.common.text_utils import get_nested as get_settings_value

try:
    from mspecviz.spectra import (
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

params = page_setup()


def file_selection_fragment():
    """File browser for idXML and mzML files."""
    _settings = st.session_state.get("settings", {})
    data_dir_default = get_settings_value(_settings, "data_directories", "default", default="/")

    idxml_dir = params.get('idxml-browser-dir', data_dir_default)
    mzml_dir = params.get('mzml-browser-dir', data_dir_default)

    st.header("File Selection")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("idXML Files (Identifications)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_idxml",
                                 extensions=['.idXML', '.idxml'],
                                 state_key='selected-idxmls',
                                 default_dir=idxml_dir, params=params)
    with col2:
        st.subheader("mzML Files (Spectra)")
        render_directory_browser(key_prefix=f"{STATE_PREFIX}_mzml",
                                 extensions=['.mzML', '.mzml'],
                                 state_key='selected-mzmls',
                                 default_dir=mzml_dir, params=params)

    for dir_param, dir_ss_key, fallback in [
        ('idxml-browser-dir', f"{STATE_PREFIX}_idxml_dir", idxml_dir),
        ('mzml-browser-dir', f"{STATE_PREFIX}_mzml_dir", mzml_dir),
    ]:
        val = st.session_state.get(dir_ss_key, fallback)
        params[dir_param] = val
        st.session_state[dir_param] = val

    save_params(params)


def reference_selection_fragment():
    """Reference substring configuration."""
    st.subheader("Reference Identification")
    mzmls = st.session_state.get('selected-mzmls', params.get('selected-mzmls', []))
    if not mzmls:
        st.warning("Please select mzML files first")
        return

    st.markdown("Enter a substring to identify reference mzML files (e.g., 'Synpep', 'Standard', 'Ref')")
    ref_substring = st.text_input(
        "Reference identifier (Optional)",
        value=params.get('reference-substring', ''),
        key="reference-substring", placeholder="e.g., Synpep"
    )

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
    available_ions = get_settings_value(st.session_state.get("settings", {}), "ion_types", "available",
                                        default=['a', 'b', 'c', 'x', 'y', 'z', 'precursor', 'immonium'])

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Ion Types**")
        current_ions = params.get('ion-types', ['b', 'y'])
        st.multiselect("Select ion types", options=available_ions,
                       default=[i for i in current_ions if i in available_ions],
                       key="ion-types")

    with col2:
        st.markdown("**m/z Tolerance**")
        st.number_input("Tolerance (Da)", min_value=0.001, max_value=1.0,
                        value=float(params.get('mz-tolerance', 0.02)),
                        step=0.005, format="%.3f", key="mz-tolerance")

    st.markdown("**Annotation Settings**")
    ann_col1, ann_col2, ann_col3 = st.columns(3)
    with ann_col1:
        st.checkbox("Add neutral losses (H2O, NH3)",
                    value=params.get('add-losses', False), key="add-losses")
    with ann_col2:
        st.checkbox("Add first prefix ion",
                    value=params.get('add-first-prefix-ion', True),
                    key="add-first-prefix-ion")
    with ann_col3:
        st.checkbox("Remove precursor peaks",
                    value=params.get('remove-precursor-peaks', False),
                    key="remove-precursor-peaks")

    st.radio("**Intensity Display Mode**",
             options=["Relative", "Raw"],
             captions=["Relative to the highest intensity", "Raw intensity"],
             index=0 if params.get('intensity-mode', 'Relative') == 'Relative' else 1,
             horizontal=True, key="intensity-mode", label_visibility="visible")

    st.markdown("**Performance**")
    st.slider("Parallel workers", min_value=1, max_value=32,
              value=params.get('workers', 4), key="workers",
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
        value=params.get("koina-enabled", koina_defaults["enabled"]),
        key="koina-enabled",
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
            current_model_id = params.get("koina-model", koina_defaults["model"])
            default_idx = model_ids.index(current_model_id) if current_model_id in model_ids else 0

            selected_model_name = st.selectbox(
                "Prediction Model", options=model_names, index=default_idx,
                key=f"{STATE_PREFIX}_koina_model_select",
                help="Deep-learning model for fragment ion prediction",
            )
            st.session_state["koina-model"] = KOINA_MODELS[selected_model_name]

        with k_col2:
            st.number_input(
                "Collision Energy (eV)", min_value=1, max_value=100,
                value=params.get("koina-ce", koina_defaults["collision_energy"]),
                key="koina-ce",
            )

        with k_col3:
            current_inst = params.get("koina-instrument", koina_defaults["instrument_type"])
            inst_idx = INSTRUMENT_TYPES.index(current_inst) if current_inst in INSTRUMENT_TYPES else 0
            st.selectbox(
                "Instrument Type", options=INSTRUMENT_TYPES, index=inst_idx,
                key="koina-instrument",
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
