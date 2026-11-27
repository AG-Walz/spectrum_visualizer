"""Explore tab — dispatcher for Spectrum Explorer / Mirror Plots / Comparison Metrics sub-views.
Registered in app.py; sub-views are @st.fragments so interactions stay scoped.
"""
import io
import logging

import matplotlib.pyplot as plt
import streamlit as st

from streamlit_app.utils.state import (
    STATE_PREFIX, get_state_key, page_setup, save_params,
)
from streamlit_app.utils.page_state import restore_analysis_if_needed, show_stale_results_warning
from core.comparison_metrics import METRIC_KEY_TO_COLUMN_NAME
from core.pipeline import COMPARISON_RT_COLUMNS
from core.plotting import (
    MEASURE_VALUE_COLUMN,
    METRIC_COLUMNS,
    SUMMARY_RENDERERS,
    peptidoform_label_series,
    prepare_single_peptide_measures,
    prepare_summary_data,
)

from streamlit_app.views.explore.explorer import render_explorer
from streamlit_app.views.explore.mirror_plots import mirror_plots_fragment
from streamlit_app.views.widgets.summary_options import render_summary_options

params = page_setup("Explore")
logger = logging.getLogger(__name__)

# Human-readable headers + one-decimal seconds for the raw RT columns of the comparison
# frame (#113). ΔRT is blank for predicted references, which carry no real retention time.
RT_COLUMN_CONFIG = {
    "query_rt": st.column_config.NumberColumn("Query RT (s)", format="%.1f"),
    "ref_rt": st.column_config.NumberColumn("Ref RT (s)", format="%.1f"),
    "rt_delta": st.column_config.NumberColumn(
        "ΔRT (s)", format="%.1f",
        help="Query RT minus reference RT; empty when the reference is predicted.",
    ),
}


# On-screen summary plots: the figures are sized for export (up to ~14 x 66 inches), so
# st.pyplot's 200 dpi made each a ~2800 x 13000 px PNG. st.image also resizes and re-encodes
# any image wider than 1460 px on every rerun, so the dpi is chosen to stay below that.
_SCREEN_MAX_WIDTH_PX = 1400
_SCREEN_MAX_DPI = 120


def _figure_png(fig) -> bytes:
    """Render a summary figure to on-screen PNG bytes and release it from pyplot."""
    dpi = min(_SCREEN_MAX_DPI, _SCREEN_MAX_WIDTH_PX / fig.get_figwidth())
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return buffer.getvalue()


# The two cached renderers below skip hashing the frame (leading underscore): results_key
# (workspace, results version) plus the filter and plot options identify it, so a rerun
# with unchanged options shows the stored picture instead of drawing it again.
@st.cache_data(show_spinner="Drawing summary plot...", max_entries=64)
def _population_plot_png(_view_frame, results_key: tuple, peptides: tuple, metric: str,
                         plot_type: str, sort_by: str) -> bytes:
    """One population plot (peptides on the y-axis) for one measure, as PNG bytes."""
    raw_frame, aggregated_frame = prepare_summary_data(_view_frame, metric=metric, sort_by=sort_by)
    return _figure_png(SUMMARY_RENDERERS[plot_type](raw_frame, aggregated_frame, metric=metric))


@st.cache_data(show_spinner="Drawing summary plot...", max_entries=64)
def _single_peptide_plot_png(_view_frame, results_key: tuple, peptide: str, metrics: tuple,
                             plot_type: str) -> bytes:
    """One plot of all selected measures (measures on the y-axis) for one peptide, as PNG
    bytes. Raises ValueError when the peptide has nothing to plot."""
    raw_frame, aggregated_frame = prepare_single_peptide_measures(_view_frame, list(metrics))
    return _figure_png(SUMMARY_RENDERERS[plot_type](
        raw_frame, aggregated_frame, metric=MEASURE_VALUE_COLUMN))


@st.fragment
def comparison_metrics_fragment():
    """Render summary plots first, then the per-PSM metrics table and per-peptidoform summary stats.
    """
    # Guard: bail when Load & Annotate hasn't produced a comparison frame
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
    if comparison_frame is None or comparison_frame.empty:
        st.info("No comparisons to show: the last run had no reference spectra. Set a reference "
                "identifier on Setup (or enable Koina), then run the analysis again.")
        return

    # One top-of-page peptide filter scoping the whole view (plots + table + stats). Empty = all.
    # Selecting exactly one peptide flips the summary plots to a measures-on-the-y-axis view (#105).
    peptidoform_labels = sorted(peptidoform_label_series(comparison_frame).unique())
    selected_peptides = st.multiselect(
        "Peptides (empty = all)", options=peptidoform_labels,
        key=f"{STATE_PREFIX}_metrics_peptide_filter",
        help="Filter every plot, table, and statistic below. Select exactly one peptide to plot all "
             "measures for it (measures on the y-axis).",
    )
    view_frame = comparison_frame
    if selected_peptides:
        view_frame = comparison_frame[peptidoform_label_series(comparison_frame).isin(selected_peptides)]
    single_peptide = selected_peptides[0] if len(selected_peptides) == 1 else None
    # Identifies the loaded results for the plot caches: cache_version holds a per-run token
    results_key = (str(st.session_state.workspace),
                   st.session_state.get(get_state_key("cache_version"), 0))

    # Summary plots: the user picks measures, plot types, and the peptide sort order via the shared
    # control (same widget the Export page uses).
    available_metrics = [m for m in METRIC_COLUMNS if m in comparison_frame.columns]
    if available_metrics:
        st.subheader("Summary Plots")
        metrics, selected_types, sort_by = render_summary_options(available_metrics)
        if not metrics or not selected_types:
            st.info("Select at least one measure and one plot type.")
        elif single_peptide:
            # One peptide: the same renderers, but with the measures on the y-axis (one row each).
            st.caption(f"All selected measures for {single_peptide}.")
            try:
                for plot_type in selected_types:
                    st.image(_single_peptide_plot_png(
                        view_frame, results_key, single_peptide, tuple(metrics), plot_type),
                        width="stretch")
            except ValueError as error:
                st.info(f"Cannot plot {single_peptide}: {error}.")
        else:
            # Several / all peptides: one population plot per measure (peptides on the y-axis).
            for metric in metrics:
                if len(metrics) > 1:
                    st.markdown(f"**{metric}**")
                for plot_type in selected_types:
                    st.image(_population_plot_png(
                        view_frame, results_key, tuple(selected_peptides), metric, plot_type, sort_by),
                        width="stretch")
        st.markdown("---")

    # Per-PSM table: one row per query/reference pair with all 7 similarity metrics
    st.subheader("Comparison Summary")
    st.caption(f"{len(view_frame)} query/reference comparisons")

    display_cols = ["Peptide", "Modification", "Charge", "Query", "Reference",
                    *COMPARISON_RT_COLUMNS,
                    "R\u00b2 match", "R\u00b2 binned", "Cosine", "Entropy",
                    "Dice", "Jaccard", "SharedRatio", "shared",
                    "query_annotation_count", "reference_annotation_count"]
    available_cols = [c for c in display_cols if c in view_frame.columns]
    st.dataframe(view_frame[available_cols], width="stretch", hide_index=True,
                 column_config=RT_COLUMN_CONFIG)

    # Summary stats: collapse PSMs into count/median/mean/std per peptidoform+charge.
    # RT columns ride along with the metrics (#113) so each peptide's elution spread and
    # its offset from the reference are visible next to the similarity scores.
    st.markdown("---")
    st.subheader("Summary Statistics")
    metric_cols = [v for v in METRIC_KEY_TO_COLUMN_NAME.values() if v in view_frame.columns]
    stat_cols = metric_cols + [c for c in COMPARISON_RT_COLUMNS if c in view_frame.columns]
    if metric_cols:
        summary = view_frame.groupby(["Peptide", "Modification", "Charge"])[stat_cols].agg(
            ["count", "median", "mean", "std"]
        )
        summary.columns = [f"{m} {s}" for m, s in summary.columns]
        st.dataframe(summary, width="stretch")


def render():
    """Render title, gate on Load & Annotate completion, then dispatch to the chosen sub-view.
    """
    # Restore analysis: rehydrate session state when the workspace has cached results
    restore_analysis_if_needed()
    st.title("Explore")
    st.caption("Interactive visualization of loaded spectra and comparison metrics")

    # Gate: Load & Annotate must have populated session state before anything renders
    if not st.session_state.get(get_state_key("data_ready")):
        st.info("No results yet. Open **Load & Annotate** and press **Run analysis**.")
        return
    show_stale_results_warning(
        params, "Press **Run analysis** on the Load & Annotate page to update them.")

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
