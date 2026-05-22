"""Rendered as the Search & Load tab — the big "Run" button plus the per-variant results table.

Heavy lifting lives in search_load_actions.py; this module handles layout, validation, and session state."""
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_param, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from core.defaults import DEFAULTS
from core.helpers import parse_peptide_list
from core.search import summarize_search_results

from streamlit_app.views.search_load.search_load_actions import (
    run_peptide_search,
    run_search_and_load,
)

params = page_setup()


def display_search_results(results: dict) -> None:
    """Render the per-variant search-results table (peptide, modification, PSM count, file resolution)."""
    summary_frame = summarize_search_results(results).collect().to_pandas()
    if not summary_frame.empty:
        display_frame = summary_frame.rename(columns={
            "base_sequence": "Peptide", "modification": "Modification",
            "psm_count": "PSMs", "idxml_files": "idXML Files",
            "mzml_resolved": "mzML Resolved", "unresolved": "Unresolved",
        })
        display_frame = display_frame.drop(columns=["display_label"], errors="ignore")
        st.dataframe(display_frame, width="stretch", hide_index=True)


@st.fragment
def search_results_fragment():
    """Render the search-results table, the Search & Load Spectra trigger, and the post-load summary."""
    # Empty state: nothing rendered until peptide search has produced variants
    variant_groups = st.session_state.get(get_state_key("variant_groups"))
    if variant_groups:
        # Section 1: per-variant results table + headline match count
        st.subheader("Search Results")
        display_search_results(variant_groups)

        total_matches = sum(len(psm_list) for psm_list in variant_groups.values())
        st.success(
            f"Found {total_matches} PSMs across {len(variant_groups)} peptide variants. "
            f"Click **Search & Load Spectra** below to load and annotate them."
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
            if st.button("Search & Load Spectra", key=f"{STATE_PREFIX}_load_btn", type="primary"):
                run_search_and_load(params)

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
        st.caption("Search results will appear here after running the search.")


def render():
    """Page entry point: rebuild analysis cache, validate inputs, then render search trigger + results fragment."""
    # Restore cached analysis (if any) and render the page header
    restore_analysis_if_needed()
    st.title("Search & Load")
    st.caption("Search for peptides in idXML files and load spectra for analysis")

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

    # Render the search trigger + results/load fragment, then persist params
    if st.button("Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(params, idxmls, peptides)

    search_results_fragment()

    save_params(params)


render()
