"""Shared summary-plot options control (measure(s), plot types, sort order).
Rendered identically on the Explore and Export pages so both offer the same choices;
the widgets share their session-state keys (summary_metrics / summary_plot_types /
summary_sort), which are the same keys core.export.RUN_SETTINGS_FIELDS reads.
"""
import streamlit as st

from core.defaults import DEFAULTS
from core.plotting import METRIC_COLUMNS
from streamlit_app.utils.page_state import SUMMARY_PLOT_LABELS, SUMMARY_SORT_LABELS


def render_summary_options(available_metrics: list[str] | None = None) -> tuple[list[str], list[str], str]:
    """Render the measure / plot-type / sort controls; return (metrics, plot_types, sort_by).
    available_metrics restricts the measure options (Explore passes only the columns present in
    the comparison frame); defaults to the full METRIC_COLUMNS list.
    """
    metric_options = list(available_metrics) if available_metrics is not None else list(METRIC_COLUMNS)
    metric_default = [DEFAULTS.summary_metric] if DEFAULTS.summary_metric in metric_options else metric_options[:1]

    # Drop any stale plot type restored from an old run (e.g. a removed "dual_axis") so the
    # multiselect's stored value never falls outside its current options.
    if "summary_plot_types" in st.session_state:
        st.session_state["summary_plot_types"] = [
            t for t in st.session_state["summary_plot_types"] if t in SUMMARY_PLOT_LABELS
        ]

    measures_col, types_col, sort_col = st.columns(3)
    with measures_col:
        metrics = st.multiselect(
            "Similarity measures", options=metric_options, default=metric_default,
            key="summary_metrics",
            help="One set of summary plots is generated per selected measure.",
        )
    with types_col:
        plot_types = st.multiselect(
            "Plot types", options=list(SUMMARY_PLOT_LABELS),
            default=list(DEFAULTS.summary_plot_types),
            format_func=lambda key: SUMMARY_PLOT_LABELS.get(key, key),
            key="summary_plot_types",
        )
    with sort_col:
        sort_by = st.selectbox(
            "Sort peptides by", options=list(SUMMARY_SORT_LABELS),
            index=list(SUMMARY_SORT_LABELS).index(DEFAULTS.summary_sort),
            format_func=lambda key: SUMMARY_SORT_LABELS.get(key, key),
            key="summary_sort",
        )
    return metrics, plot_types, sort_by
