"""Rendered as the Load & Annotate tab — the big "Run" button plus the per-peptidoform results table.

Heavy lifting lives in load_annotate_actions.py; this module handles layout, validation, and session state."""
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_param, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from core.defaults import DEFAULTS
from core.helpers import parse_peptide_list
from core.peptide_lookup import summarize_lookup_results

from streamlit_app.views.load_annotate.load_annotate_actions import (
    run_peptide_lookup,
    run_load_and_annotate,
)

params = page_setup("Load & Annotate")


def display_lookup_results(peptidoform_psms: dict) -> None:
    """Render the per-peptidoform lookup-results table (peptide, modification, PSM count, file resolution)."""
    summary_frame = summarize_lookup_results(peptidoform_psms).collect().to_pandas()
    if not summary_frame.empty:
        display_frame = summary_frame.rename(columns={
            "base_sequence": "Peptide", "modification": "Modification",
            "psm_count": "PSMs", "idxml_files": "idXML Files",
            "mzml_resolved": "mzML Resolved", "unresolved": "Unresolved",
        })
        display_frame = display_frame.drop(columns=["display_label"], errors="ignore")
        st.dataframe(display_frame, width="stretch", hide_index=True)


@st.fragment
def lookup_results_fragment():
    """Render the lookup-results table, the Load & Annotate Spectra trigger, and the post-load summary."""
    # Empty state: nothing rendered until peptide lookup has produced peptidoforms
    peptidoform_psms = st.session_state.get(get_state_key("peptidoform_psms"))
    if peptidoform_psms:
        # Section 1: per-peptidoform results table + headline match count
        st.subheader("Lookup Results")
        display_lookup_results(peptidoform_psms)

        total_matches = sum(len(psm_list) for psm_list in peptidoform_psms.values())
        st.success(
            f"Found {total_matches} PSMs across {len(peptidoform_psms)} peptidoforms. "
            f"Click **Load & Annotate Spectra** below to load and annotate them."
        )
        # Section 2: load-spectra controls — only when mzMLs were picked on the Setup page
        mzmls = get_param(params, "selected_mzmls", [])
        if mzmls:
            st.subheader("Load Spectra & Compute Metrics")
            # Hint: what the loader will do given the current reference / Koina settings
            reference_substring = get_param(params, "reference_substring", "")
            if reference_substring:
                st.caption(
                    f"Will pair query spectra with references (substring: '{reference_substring}') "
                    f"and compute 7 comparison metrics per pair."
                )
            else:
                koina_enabled = get_param(params, "koina_enabled", DEFAULTS.koina_enabled)
                if koina_enabled:
                    st.caption("No reference configured -- will use Koina AI prediction as fallback.")
                else:
                    st.caption("No reference configured -- will annotate individual spectra.")
            if st.button("Load & Annotate Spectra", key=f"{STATE_PREFIX}_load_btn", type="primary"):
                run_load_and_annotate(params)

            # Section 3: post-load summary — only after the action has set data_ready
            if st.session_state.get(get_state_key("data_ready")):
                comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
                if comparison_frame is not None and not comparison_frame.empty:
                    st.success(
                        f"Data loaded with {len(comparison_frame)} query/reference comparisons. "
                        f"Go to **Explore** page."
                    )
                else:
                    st.success("Data loaded. Go to **Explore** page.")
        else:
            st.warning("Please select mzML files first (Setup page)")
    else:
        st.caption("Lookup results will appear here after running the lookup.")


def render():
    """Page entry point: rebuild analysis cache, validate inputs, then render lookup trigger + results fragment."""
    # Restore cached analysis (if any) and render the page header
    restore_analysis_if_needed()
    st.title("Load & Annotate")
    st.caption("Look up peptides in idXML files and load spectra for analysis")

    # Resolve inputs from session params and bail early if Setup hasn't been completed
    idxmls = get_param(params, "selected_idxmls", [])
    peptide_text = get_param(params, "target_peptides", "")
    peptides = parse_peptide_list(peptide_text)

    if not idxmls:
        st.warning("Please select idXML files first (Setup page)")
        return
    if not peptides:
        st.warning("Please enter target peptides first (Setup page)")
        return

    # Render the lookup trigger + results/load fragment, then persist params
    if st.button("Look Up Peptides", key=f"{STATE_PREFIX}_lookup_btn", type="primary"):
        run_peptide_lookup(params, idxmls, peptides)

    lookup_results_fragment()

    save_params(params)


render()
