"""Explore tab — dispatcher for Spectrum Explorer / Mirror Plots / Comparison Metrics sub-views.
Registered in app.py; sub-views are @st.fragments so interactions stay scoped.
"""
import logging

import matplotlib.pyplot as plt
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed
from core.comparison_metrics import METRIC_KEY_TO_COLUMN_NAME
from core.defaults import DEFAULTS
from core.plotting import (
    METRIC_COLUMNS,
    plot_classic_dual_axis,
    plot_ridgeline,
    prepare_summary_data,
)

from streamlit_app.views.explore.explorer import render_explorer
from streamlit_app.views.explore.mirror_plots import mirror_plots_fragment

params = page_setup()
logger = logging.getLogger(__name__)


@st.fragment
def comparison_metrics_fragment():
    """Render per-PSM metrics table, per-peptidoform summary stats, and summary plots.
    """
    # Guard: bail when Load & Annotate hasn't produced a comparison frame
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
    if comparison_frame is None or comparison_frame.empty:
        st.info("No comparison metrics available. Configure a reference substring and re-run.")
        return

    # Per-PSM table: one row per query/reference pair with all 7 similarity metrics
    st.subheader("Comparison Summary")
    st.caption(f"{len(comparison_frame)} query/reference comparisons")

    display_cols = ["Peptide", "Modification", "Charge", "Query", "Reference",
                    "R\u00b2 match", "R\u00b2 binned", "Cosine", "Entropy",
                    "Dice", "Jaccard", "SharedRatio", "shared",
                    "query_annotation_count", "reference_annotation_count"]
    available_cols = [c for c in display_cols if c in comparison_frame.columns]
    st.dataframe(comparison_frame[available_cols], width="stretch", hide_index=True)

    # Summary stats: collapse PSMs into count/median/mean/std per peptidoform+charge
    st.markdown("---")
    st.subheader("Summary Statistics")
    metric_cols = [v for v in METRIC_KEY_TO_COLUMN_NAME.values() if v in comparison_frame.columns]
    if metric_cols:
        summary = comparison_frame.groupby(["Peptide", "Modification", "Charge"])[metric_cols].agg(
            ["count", "median", "mean", "std"]
        )
        summary.columns = [f"{m} {s}" for m, s in summary.columns]
        st.dataframe(summary, width="stretch")

    # Summary plots: dual-axis bar/line + ridgeline density for the chosen metric
    available_metrics = [m for m in METRIC_COLUMNS if m in comparison_frame.columns]
    if available_metrics:
        st.markdown("---")
        st.subheader("Summary Plots")
        metric = st.selectbox(
            "Metric", options=available_metrics,
            index=available_metrics.index(DEFAULTS.summary_metric) if DEFAULTS.summary_metric in available_metrics else 0,
            key=f"{STATE_PREFIX}_summary_metric",
        )
        raw_frame, aggregated_frame = prepare_summary_data(comparison_frame, metric=metric)
        fig_dual = plot_classic_dual_axis(raw_frame, aggregated_frame, metric=metric)
        st.pyplot(fig_dual)
        plt.close(fig_dual)
        fig_ridge = plot_ridgeline(raw_frame, aggregated_frame, metric=metric)
        st.pyplot(fig_ridge)
        plt.close(fig_ridge)


def render():
    """Render title, gate on Load & Annotate completion, then dispatch to the chosen sub-view.
    """
    # Restore analysis: rehydrate session state when the workspace has cached results
    restore_analysis_if_needed()
    st.title("Explore")
    st.caption("Interactive visualization of loaded spectra and comparison metrics")

    # Gate: Load & Annotate must have populated session state before anything renders
    if not st.session_state.get(get_state_key("data_ready")):
        st.info("Run 'Load & Annotate' first to populate the explorer.")
        return

    # Sub-view picker: horizontal radio that survives reruns via STATE_PREFIX key
    explore_mode = st.radio(
        "View mode", options=["Spectrum Explorer", "Mirror Plots", "Comparison Metrics"],
        horizontal=True, key=f"{STATE_PREFIX}_explore_mode"
    )

    # Dispatch: only the selected sub-view's fragment runs

    if explore_mode == "Spectrum Explorer":
        render_explorer(params)
    elif explore_mode == "Mirror Plots":
        mirror_plots_fragment(params)
    elif explore_mode == "Comparison Metrics":
        comparison_metrics_fragment()

    save_params(params)


render()
