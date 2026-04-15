"""
Search & Load Page — Streamlit entry point.

Orchestrates the page UI (buttons, fragments, results table). Long-running
ops live in ``pipeline.py``; this file stays focused on layout and state.
"""
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_param, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from mspecviz.common.constants import DEFAULTS
from mspecviz.common.text_utils import parse_peptide_list
from mspecviz.search import summarize_search_results

from streamlit_app.views.search_load.pipeline import (
    run_peptide_search,
    run_search_and_load,
)

params = page_setup()


def display_search_results(results: dict) -> None:
    """Display search results table."""
    summary_lf = summarize_search_results(results)
    summary_df = summary_lf.collect().to_pandas()
    if not summary_df.empty:
        display_df = summary_df.rename(columns={
            "base_sequence": "Peptide", "modification": "Modification",
            "psm_count": "PSMs", "idxml_files": "idXML Files",
            "mzml_resolved": "mzML Resolved", "unresolved": "Unresolved",
        })
        display_df = display_df.drop(columns=["display_label"], errors="ignore")
        st.dataframe(display_df, width="stretch", hide_index=True)


@st.fragment
def search_results_fragment():
    """Display search results and load controls as a fragment."""
    variant_groups = st.session_state.get(get_state_key("variant_groups"))
    if variant_groups:
        st.subheader("Search Results")
        display_search_results(variant_groups)

        total_matches = sum(len(variant_group.matches) for variant_group in variant_groups.values())
        st.success(
            f"Found {total_matches} PSMs across {len(variant_groups)} peptide variants. "
            f"Click **Search & Load Spectra** below to load and annotate them."
        )
        mzmls = get_param(params, "selected-mzmls", DEFAULTS["selected-mzmls"])
        if mzmls:
            st.subheader("Load Spectra & Compute Metrics")
            reference_substring = get_param(params, "reference-substring", DEFAULTS["reference-substring"])
            if reference_substring:
                st.caption(
                    f"Will pair query spectra with references (substring: '{reference_substring}') "
                    f"and compute 7 comparison metrics per pair."
                )
            else:
                koina_enabled = get_param(params, "koina-enabled", DEFAULTS["koina-enabled"])
                if koina_enabled:
                    st.caption("No reference configured -- will use Koina AI prediction as fallback.")
                else:
                    st.caption("No reference configured -- will annotate individual spectra.")
            if st.button("Search & Load Spectra", key=f"{STATE_PREFIX}_load_btn", type="primary"):
                run_search_and_load(params)

            if st.session_state.get(get_state_key("data_ready")):
                comparison_df = st.session_state.get(get_state_key("comparison_df"))
                if comparison_df is not None and not comparison_df.empty:
                    st.success(
                        f"Data loaded with {len(comparison_df)} query/reference comparisons. "
                        f"Go to **Explore** page."
                    )
                else:
                    st.success("Data loaded. Go to **Explore** page.")
        else:
            st.warning("Please select mzML files first (Setup page)")
    else:
        st.caption("Search results will appear here after running the search.")


def render():
    restore_analysis_if_needed()
    st.title("Search & Load")
    st.caption("Search for peptides in idXML files and load spectra for analysis")

    idxmls = get_param(params, "selected-idxmls", DEFAULTS["selected-idxmls"])
    peptide_text = get_param(params, "target-peptides", DEFAULTS["target-peptides"])
    peptides = parse_peptide_list(peptide_text)

    if not idxmls:
        st.warning("Please select idXML files first (Setup page)")
        return
    if not peptides:
        st.warning("Please enter target peptides first (Setup page)")
        return

    if st.button("Search Peptides", key=f"{STATE_PREFIX}_search_btn", type="primary"):
        run_peptide_search(params, idxmls, peptides)

    search_results_fragment()

    save_params(params)


render()
