"""
Explore Page — Streamlit entry point with 3 sub-views.

Each sub-view is an @st.fragment so interactions within one view
(e.g., selecting a PSM in mirror plots) don't rerun the others.
"""
import logging

import matplotlib.pyplot as plt
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from mspecviz.annotation import METRIC_SHORT_TO_DISPLAY
from mspecviz.plotting import (
    DEFAULT_METRIC,
    METRIC_COLUMNS,
    plot_classic_dual_axis,
    plot_ridgeline,
    prepare_summary_data,
)

from streamlit_app.views.explore.explorer import openms_insight_explorer_fragment
from streamlit_app.views.explore.plots import mirror_plots_fragment

params = page_setup()
logger = logging.getLogger(__name__)


@st.fragment
def comparison_metrics_fragment():
    """Show comparison metrics table and summary statistics."""
    comparison_df = st.session_state.get(get_state_key("comparison_df"))
    if comparison_df is None or comparison_df.empty:
        st.info("No comparison metrics available. Configure a reference substring and re-run.")
        return

    st.subheader("Comparison Summary")
    st.caption(f"{len(comparison_df)} query/reference comparisons")

    display_cols = ["Peptide", "Modification", "Charge", "Query", "Reference",
                    "R\u00b2 match", "R\u00b2 binned", "Cosine", "Entropy",
                    "Dice", "Jaccard", "SharedRatio", "shared", "ann_q", "ann_r"]
    available_cols = [c for c in display_cols if c in comparison_df.columns]
    st.dataframe(comparison_df[available_cols], width="stretch", hide_index=True)

    # Summary statistics per peptidoform
    st.markdown("---")
    st.subheader("Summary Statistics")
    import pandas as pd  # noqa: F401 — reserved for future formatting tweaks
    metric_cols = [v for v in METRIC_SHORT_TO_DISPLAY.values() if v in comparison_df.columns]
    if metric_cols:
        summary = comparison_df.groupby(["Peptide", "Modification", "Charge"])[metric_cols].agg(
            ["count", "median", "mean", "std"]
        )
        summary.columns = [f"{m} {s}" for m, s in summary.columns]
        st.dataframe(summary, width="stretch")

    # Summary plots (same as export)
    available_metrics = [m for m in METRIC_COLUMNS if m in comparison_df.columns]
    if available_metrics:
        st.markdown("---")
        st.subheader("Summary Plots")
        metric = st.selectbox(
            "Metric", options=available_metrics,
            index=available_metrics.index(DEFAULT_METRIC) if DEFAULT_METRIC in available_metrics else 0,
            key=f"{STATE_PREFIX}_summary_metric",
        )
        try:
            raw_df, agg_df = prepare_summary_data(comparison_df, metric=metric)
            fig_dual = plot_classic_dual_axis(raw_df, agg_df, metric=metric)
            st.pyplot(fig_dual)
            plt.close(fig_dual)
            fig_ridge = plot_ridgeline(raw_df, agg_df, metric=metric)
            st.pyplot(fig_ridge)
            plt.close(fig_ridge)
        except Exception as e:
            st.warning(f"Could not generate summary plots: {e}")


def render():
    restore_analysis_if_needed()
    st.title("Explore")
    st.caption("Interactive visualization of loaded spectra and comparison metrics")

    if not st.session_state.get(get_state_key("data_ready")):
        st.info("Run 'Search & Load' first to populate the explorer.")
        return

    explore_mode = st.radio(
        "View mode", options=["Interactive Explorer", "Mirror Plots", "Comparison Metrics"],
        horizontal=True, key=f"{STATE_PREFIX}_explore_mode"
    )

    if explore_mode == "Interactive Explorer":
        openms_insight_explorer_fragment(params)
    elif explore_mode == "Mirror Plots":
        mirror_plots_fragment(params)
    elif explore_mode == "Comparison Metrics":
        comparison_metrics_fragment()

    save_params(params)


render()
