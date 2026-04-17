"""
Interactive Plotly mirror-plot fragment.

Lets the user pick any query PSM from a dropdown; if a reference is available
in ``reference_store`` for its variant, a mirror plot (+ metrics table) is
rendered, otherwise a single-spectrum view.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from streamlit_app.utils.state import STATE_PREFIX, get_param, get_state_key
from mspecviz.annotation import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    format_metrics_text,
)
from mspecviz.common.defaults import DEFAULTS
from mspecviz.pipeline import build_psm_to_variant_mapping, format_spectrum_labels
from mspecviz.plotting import create_plotly_mirror_plot, create_plotly_single_plot


@st.fragment
def mirror_plots_fragment(params: dict):
    """Interactive Plotly mirror plots for each PSM vs its references."""
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    reference_store = st.session_state.get(get_state_key("reference_store"), {})
    variant_groups = st.session_state.get(get_state_key("variant_groups"), {})

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

    selected_label = st.selectbox("Select PSM", options=[label for _, label in psm_options],
                                  key=f"{STATE_PREFIX}_mirror_psm_select")

    selected_psm_id = None
    for psm_id, label in psm_options:
        if label == selected_label:
            selected_psm_id = psm_id
            break

    if selected_psm_id is None or selected_psm_id not in annotated_queries:
        return

    query_spectrum, query_frame = annotated_queries[selected_psm_id]
    variant_key = psm_id_to_variant_key.get(selected_psm_id)
    references = reference_store.get(variant_key, []) if variant_key else []

    normalize = get_param(params, "intensity_mode", DEFAULTS.intensity_mode) == "Relative"

    if references:
        for reference_spectrum, reference_frame in references:
            metrics = compute_all_metrics(query_frame, reference_frame)
            metrics_text = format_metrics_text(metrics)

            query_label, ref_label = format_spectrum_labels(query_spectrum, reference_spectrum)
            title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"

            fig = create_plotly_mirror_plot(
                upper_spectrum=query_frame,
                lower_spectrum=reference_frame,
                title=title,
                upper_label=query_label,
                lower_label=ref_label,
                metrics_text=metrics_text,
                normalize=normalize,
            )
            st.plotly_chart(fig, width="stretch")

            with st.expander("Comparison Metrics", expanded=False):
                metric_data = {
                    METRIC_SHORT_TO_DISPLAY[k]: f"{v:.4f}" if not (isinstance(v, float) and v != v) else "N/A"
                    for k, v in metrics.items() if k in METRIC_SHORT_TO_DISPLAY
                }
                metric_data["Shared Ions"] = str(metrics.get("shared", 0))
                metric_data["Annotated (Query)"] = str(metrics.get("ann_q", 0))
                metric_data["Annotated (Ref)"] = str(metrics.get("ann_r", 0))
                st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)
    else:
        title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"
        label = f"{Path(query_spectrum.source_file).stem} (RT: {query_spectrum.rt:.1f})"
        fig = create_plotly_single_plot(
            spectrum=query_frame, title=title, label=label, normalize=normalize,
        )
        st.plotly_chart(fig, width="stretch")
