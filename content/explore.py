"""
Explore Page — Interactive openms_insight components, Plotly mirror plots, comparison metrics.

Each sub-view is an @st.fragment so interactions within one view
(e.g., selecting a PSM in mirror plots) don't rerun the others.
"""
import os
import logging
import streamlit as st
import polars as pl
from pathlib import Path

from src.common.common import page_setup, save_params
from content._shared import (
    STATE_PREFIX, get_state_key, restore_analysis_if_needed,
)

params = page_setup()
logger = logging.getLogger(__name__)


def _p(key, default=None):
    """Get a param value: session_state first (if widget present), then params.json."""
    return st.session_state.get(key, params.get(key, default))


# ---------------------------------------------------------------------------
# Sub-view: openms_insight Interactive Explorer
# ---------------------------------------------------------------------------

@st.fragment
def openms_insight_explorer_fragment():
    """openms_insight Table + LinePlot + SequenceView + Heatmap."""
    from openms_insight import Table, LinePlot, SequenceView, Heatmap, StateManager

    _export_dir = _p('export-dir', '')
    CACHE_PATH = _export_dir if _export_dir else str(Path(st.session_state.workspace) / "exports")
    os.makedirs(CACHE_PATH, exist_ok=True)

    cache_ver = st.session_state.get(get_state_key('cache_version'), 0)
    regen_key = get_state_key('cache_rendered_version')
    needs_regen = st.session_state.get(regen_key, -1) != cache_ver

    if 'msinsight_sm' not in st.session_state:
        st.session_state['msinsight_sm'] = StateManager(session_key="msinsight_state")
    sm = st.session_state['msinsight_sm']

    psm_lf = st.session_state.get(get_state_key('psm_lf'))
    peaks_lf = st.session_state.get(get_state_key('peaks_lf'))
    seq_lf = st.session_state.get(get_state_key('seq_lf'))
    heatmap_lf = st.session_state.get(get_state_key('heatmap_lf'))

    # Show currently selected PSM info
    current_state = sm.get_state_for_vue()
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
    if needs_regen:
        table = Table(
            cache_id="psm_table",
            data=psm_lf,
            interactivity={"selected_psm": "id"},
            column_definitions=psm_columns,
            cache_path=CACHE_PATH,
            regenerate_cache=True,
            index_field="id",
        )
    else:
        table = Table(cache_id="psm_table", cache_path=CACHE_PATH)
    table(key="psm_table", state_manager=sm)

    # SequenceView
    st.subheader("Sequence Coverage")
    if needs_regen:
        peaks_for_sv = peaks_lf.with_row_index("peak_id").rename({"x": "mass", "y": "intensity"})
        sv = SequenceView(
            cache_id="seq_view",
            sequence_data=seq_lf,
            peaks_data=peaks_for_sv,
            filters={"selected_psm": "psm_id"},
            annotation_config={
                "ion_types": _p('ion-types', ['b', 'y']),
                "tolerance": _p('mz-tolerance', 0.02),
                "tolerance_ppm": False,
            },
            cache_path=CACHE_PATH,
            height=350,
        )
    else:
        sv = SequenceView(cache_id="seq_view", cache_path=CACHE_PATH)
    sv(key="seq_view", state_manager=sm)

    # LinePlot linked to SequenceView
    st.subheader("Fragment Spectrum")
    if needs_regen:
        lineplot = LinePlot(
            cache_id="spectrum_plot",
            data=peaks_lf,
            filters={"selected_psm": "psm_id"},
            filter_defaults={"selected_psm": 0},
            cache_path=CACHE_PATH,
            regenerate_cache=True,
            x_column="x", y_column="y",
            highlight_column="highlight",
            annotation_column="annotation",
            x_label="m/z", y_label="Intensity",
        )
    else:
        lineplot = LinePlot(cache_id="spectrum_plot", cache_path=CACHE_PATH)
    lineplot(key="spectrum_plot", state_manager=sm, sequence_view_key="seq_view")

    # Heatmap
    st.subheader("RT x m/z Heatmap")
    if needs_regen:
        heatmap_comp = Heatmap(
            cache_id="rt_heatmap",
            data=heatmap_lf,
            x_column="x", y_column="y", intensity_column="intensity",
            x_label="Retention Time (s)", y_label="m/z",
            log_scale=True,
            cache_path=CACHE_PATH,
            regenerate_cache=True,
        )
    else:
        heatmap_comp = Heatmap(cache_id="rt_heatmap", cache_path=CACHE_PATH)
    heatmap_comp(key="rt_heatmap", state_manager=sm)

    if needs_regen:
        st.session_state[regen_key] = cache_ver


# ---------------------------------------------------------------------------
# Sub-view: Mirror Plots
# ---------------------------------------------------------------------------

@st.fragment
def mirror_plots_fragment():
    """Interactive Plotly mirror plots for each PSM vs its references."""
    from src.mirror_plot import (
        create_plotly_mirror_plot, create_plotly_single_plot, ION_COLORS,
    )
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text

    annotated_queries = st.session_state.get(get_state_key('annotated_queries'), {})
    reference_store = st.session_state.get(get_state_key('reference_store'), {})
    variant_groups = st.session_state.get(get_state_key('variant_groups'), {})

    if not annotated_queries:
        st.info("No annotated spectra available.")
        return

    # Build list of PSMs for selection
    psm_options = []
    psm_id_to_vk = {}
    for vk, vg in variant_groups.items():
        for match in vg.matches:
            from src.spectrum_loader import spectrum_lookup_key
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for pid, (spec, _) in annotated_queries.items():
                if spectrum_lookup_key(spec) == key:
                    label = f"[{pid}] {vg.display_label} -- {Path(spec.source_file).stem} (RT: {spec.rt:.1f})"
                    psm_options.append((pid, label))
                    psm_id_to_vk[pid] = vk
                    break

    if not psm_options:
        for pid, (spec, _) in sorted(annotated_queries.items()):
            label = f"[{pid}] {spec.peptidoform} -- {Path(spec.source_file).stem}"
            psm_options.append((pid, label))

    selected_label = st.selectbox("Select PSM", options=[lbl for _, lbl in psm_options],
                                  key=f"{STATE_PREFIX}_mirror_psm_select")

    selected_pid = None
    for pid, lbl in psm_options:
        if lbl == selected_label:
            selected_pid = pid
            break

    if selected_pid is None or selected_pid not in annotated_queries:
        return

    query_spec, query_df = annotated_queries[selected_pid]
    vk = psm_id_to_vk.get(selected_pid)
    refs = reference_store.get(vk, []) if vk else []

    normalize = _p('intensity-mode', 'Relative') == 'Relative'

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
                import pandas as pd
                from src.spectrum_metrics import METRIC_SHORT_TO_DISPLAY
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


# ---------------------------------------------------------------------------
# Sub-view: Comparison Metrics
# ---------------------------------------------------------------------------

@st.fragment
def comparison_metrics_fragment():
    """Show comparison metrics table and summary statistics."""
    import pandas as pd

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

    # Ion match statistics
    ims_cols = [c for c in comparison_df.columns if c.startswith("ims_")]
    if ims_cols:
        with st.expander(f"Ion Match Statistics ({len(ims_cols)} metrics)", expanded=False):
            ims_display_cols = ["Peptide", "Query", "Reference"] + ims_cols
            ims_available = [c for c in ims_display_cols if c in comparison_df.columns]
            st.dataframe(comparison_df[ims_available], width="stretch", hide_index=True)

    # Summary statistics per peptidoform
    st.markdown("---")
    st.subheader("Summary Statistics")
    from src.spectrum_metrics import METRIC_SHORT_TO_DISPLAY
    metric_cols = [v for v in METRIC_SHORT_TO_DISPLAY.values() if v in comparison_df.columns]
    if metric_cols:
        summary = comparison_df.groupby(["Peptide", "Modification", "Charge"])[metric_cols].agg(
            ["count", "median", "mean", "std"]
        )
        summary.columns = [f"{m} {s}" for m, s in summary.columns]
        st.dataframe(summary, width="stretch")


# ---------------------------------------------------------------------------
# Page render
# ---------------------------------------------------------------------------

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
