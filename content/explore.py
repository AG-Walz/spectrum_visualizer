"""
Explore Page — Interactive openms_insight components, Plotly mirror plots, comparison metrics.

Each sub-view is an @st.fragment so interactions within one view
(e.g., selecting a PSM in mirror plots) don't rerun the others.
"""
import os
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import polars as pl
from openms_insight import Table, LinePlot, SequenceView, Heatmap, StateManager

from src.common.streamlit_setup import page_setup, save_params
from src.mirror_plot import (
    create_plotly_mirror_plot, create_plotly_single_plot, ION_COLORS,
)
from src.spectrum_metrics import compute_all_metrics, format_metrics_text, METRIC_SHORT_TO_DISPLAY
from src.summary_plots import (
    METRIC_COLUMNS, DEFAULT_METRIC,
    prepare_summary_data, plot_classic_dual_axis, plot_ridgeline,
)
from src.query_processing import build_psm_to_variant_mapping
from content.page_state import (
    STATE_PREFIX, get_state_key, get_param, restore_analysis_if_needed,
)

params = page_setup()
logger = logging.getLogger(__name__)


def _create_component(cls, cache_id, cache_path, needs_regeneration, **kwargs):
    """Create an openms_insight component, using cache when available."""
    if needs_regeneration:
        return cls(cache_id=cache_id, cache_path=cache_path, **kwargs)
    return cls(cache_id=cache_id, cache_path=cache_path)


@st.fragment
def openms_insight_explorer_fragment():
    """openms_insight Table + LinePlot + SequenceView + Heatmap."""
    export_dir = get_param(params, 'export-dir', '')
    CACHE_PATH = export_dir if export_dir else str(Path(st.session_state.workspace) / "exports")
    os.makedirs(CACHE_PATH, exist_ok=True)

    cache_version = st.session_state.get(get_state_key('cache_version'), 0)
    regeneration_key = get_state_key('cache_rendered_version')
    needs_regeneration = st.session_state.get(regeneration_key, -1) != cache_version

    if 'msinsight_sm' not in st.session_state:
        st.session_state['msinsight_sm'] = StateManager(session_key="msinsight_state")
    state_manager = st.session_state['msinsight_sm']

    psm_lf = st.session_state.get(get_state_key('psm_lf'))
    peaks_lf = st.session_state.get(get_state_key('peaks_lf'))
    seq_lf = st.session_state.get(get_state_key('seq_lf'))
    heatmap_lf = st.session_state.get(get_state_key('heatmap_lf'))

    # Show currently selected PSM info
    current_state = state_manager.get_state_for_vue()
    selected_id = current_state.get("selected_psm")
    if selected_id is not None and psm_lf is not None:
        try:
            row = psm_lf.filter(pl.col("id") == selected_id).collect()
            if row.height > 0:
                label = row["display_label"][0]
                source = row["source_file"][0]
                st.info(f"Viewing: **{label}** (PSM {selected_id}) from {source}")
        except Exception:
            pass

    # PSM Table
    st.subheader("PSM Table")
    st.caption("Click a row to update the spectrum and sequence views below.")
    psm_columns = [
        {"field": "display_label", "title": "Peptide"},
        {"field": "charge", "title": "Z"},
        {"field": "score", "title": "Score", "formatter": "money", "formatterParams": {"precision": 4}},
        {"field": "source_file", "title": "Source File"},
    ]
    table = _create_component(Table, "psm_table", CACHE_PATH, needs_regeneration,
                              data=psm_lf, interactivity={"selected_psm": "id"},
                              column_definitions=psm_columns, regenerate_cache=True, index_field="id")
    table(key="psm_table", state_manager=state_manager)

    # SequenceView
    st.subheader("Sequence Coverage")
    if needs_regeneration:
        peaks_for_sv = peaks_lf.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
    else:
        peaks_for_sv = None
    sequence_view = _create_component(
        SequenceView, "seq_view", CACHE_PATH, needs_regeneration,
        sequence_data=seq_lf, peaks_data=peaks_for_sv,
        filters={"selected_psm": "psm_id"},
        annotation_config={
            "ion_types": get_param(params, 'ion-types', ['b', 'y']),
            "tolerance": get_param(params, 'mz-tolerance', 0.02),
            "tolerance_ppm": False,
        },
        height=350,
    )
    sequence_view(key="seq_view", state_manager=state_manager)

    # LinePlot linked to SequenceView
    st.subheader("Fragment Spectrum")
    lineplot = _create_component(
        LinePlot, "spectrum_plot", CACHE_PATH, needs_regeneration,
        data=peaks_lf, filters={"selected_psm": "psm_id"},
        filter_defaults={"selected_psm": 0}, regenerate_cache=True,
        x_column="x", y_column="y", highlight_column="highlight",
        annotation_column="annotation", x_label="m/z", y_label="Intensity",
    )
    lineplot(key="spectrum_plot", state_manager=state_manager, sequence_view_key="seq_view")

    # Heatmap
    st.subheader("RT x m/z Heatmap")
    heatmap_comp = _create_component(
        Heatmap, "rt_heatmap", CACHE_PATH, needs_regeneration,
        data=heatmap_lf, x_column="x", y_column="y", intensity_column="intensity",
        x_label="Retention Time (s)", y_label="m/z", log_scale=True, regenerate_cache=True,
    )
    heatmap_comp(key="rt_heatmap", state_manager=state_manager)

    if needs_regeneration:
        st.session_state[regeneration_key] = cache_version


@st.fragment
def mirror_plots_fragment():
    """Interactive Plotly mirror plots for each PSM vs its references."""
    annotated_queries = st.session_state.get(get_state_key('annotated_queries'), {})
    reference_store = st.session_state.get(get_state_key('reference_store'), {})
    variant_groups = st.session_state.get(get_state_key('variant_groups'), {})

    if not annotated_queries:
        st.info("No annotated spectra available.")
        return

    # Build list of PSMs for selection
    psm_id_to_variant_key = build_psm_to_variant_mapping(variant_groups, annotated_queries)
    psm_options = []
    for psm_id, variant_key in psm_id_to_variant_key.items():
        spec, _ = annotated_queries[psm_id]
        variant_group = variant_groups[variant_key]
        label = f"[{psm_id}] {variant_group.display_label} -- {Path(spec.source_file).stem} (RT: {spec.rt:.1f})"
        psm_options.append((psm_id, label))

    if not psm_options:
        for psm_id, (spec, _) in sorted(annotated_queries.items()):
            label = f"[{psm_id}] {spec.peptidoform} -- {Path(spec.source_file).stem}"
            psm_options.append((psm_id, label))

    selected_label = st.selectbox("Select PSM", options=[lbl for _, lbl in psm_options],
                                  key=f"{STATE_PREFIX}_mirror_psm_select")

    selected_psm_id = None
    for psm_id, lbl in psm_options:
        if lbl == selected_label:
            selected_psm_id = psm_id
            break

    if selected_psm_id is None or selected_psm_id not in annotated_queries:
        return

    query_spec, query_df = annotated_queries[selected_psm_id]
    variant_key = psm_id_to_variant_key.get(selected_psm_id)
    refs = reference_store.get(variant_key, []) if variant_key else []

    normalize = get_param(params, 'intensity-mode', 'Relative') == 'Relative'

    if refs:
        for ref_spec, ref_df in refs:
            metrics = compute_all_metrics(query_df, ref_df)
            metrics_text = format_metrics_text(metrics)

            query_label = f"Query: {Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
            ref_source = Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina"
            ref_label = f"Ref: {ref_source} (RT: {ref_spec.rt:.1f})"
            title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"

            fig = create_plotly_mirror_plot(
                upper_spectrum=query_df,
                lower_spectrum=ref_df,
                title=title,
                upper_label=query_label,
                lower_label=ref_label,
                metrics_text=metrics_text,
                normalize=normalize,
            )
            st.plotly_chart(fig, width="stretch")

            with st.expander("Comparison Metrics", expanded=False):
                metric_data = {METRIC_SHORT_TO_DISPLAY[k]: f"{v:.4f}" if not (isinstance(v, float) and v != v) else "N/A"
                               for k, v in metrics.items() if k in METRIC_SHORT_TO_DISPLAY}
                metric_data["Shared Ions"] = str(metrics.get("shared", 0))
                metric_data["Annotated (Query)"] = str(metrics.get("ann_q", 0))
                metric_data["Annotated (Ref)"] = str(metrics.get("ann_r", 0))
                st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)
    else:
        title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
        label = f"{Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
        fig = create_plotly_single_plot(
            spectrum=query_df, title=title, label=label, normalize=normalize,
        )
        st.plotly_chart(fig, width="stretch")


@st.fragment
def comparison_metrics_fragment():
    """Show comparison metrics table and summary statistics."""
    comparison_df = st.session_state.get(get_state_key('comparison_df'))
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

    if not st.session_state.get(get_state_key('data_ready')):
        st.info("Run 'Search & Load' first to populate the explorer.")
        return

    explore_mode = st.radio(
        "View mode", options=["Interactive Explorer", "Mirror Plots", "Comparison Metrics"],
        horizontal=True, key=f"{STATE_PREFIX}_explore_mode"
    )

    if explore_mode == "Interactive Explorer":
        openms_insight_explorer_fragment()
    elif explore_mode == "Mirror Plots":
        mirror_plots_fragment()
    elif explore_mode == "Comparison Metrics":
        comparison_metrics_fragment()

    save_params(params)


render()
