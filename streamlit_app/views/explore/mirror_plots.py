"""Mirror Plots sub-view of the Explore tab — Plotly mirror plot driven by a two-stage picker.
Rendered by explore_page.py: peptide selectbox + two Tables (query left, reference right).
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import StateManager, Table

from streamlit_app.utils.state import STATE_PREFIX, get_param, get_state_key
from core.annotation import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    format_metrics_text,
)
from core.defaults import DEFAULTS
from core.pipeline import (
    build_psm_to_peptidoform_mapping,
    format_spectrum_labels,
)
from core.plotting import create_spectrum_plot
from core.spectra.streaming_types import source_label
from core.spectra.mzml_native_ids import extract_id_number


def _eligible_peptidoforms(peptidoform_psms: dict, reference_store: dict, psm_id_to_peptidoform: dict) -> list:
    """Return peptidoforms that have at least one query PSM AND at least one reference loaded.
    """
    peptidoforms_with_query = set(psm_id_to_peptidoform.values())
    return [peptidoform for peptidoform in peptidoform_psms
            if peptidoform in peptidoforms_with_query and reference_store.get(peptidoform)]


def _ref_row_id(peptidoform: str, position: int) -> int:
    """Hash (peptidoform, position) into a unique row id so the reference Table refreshes per peptide.
    Plain positional ids collide across peptides; salting with peptidoform busts the Vue render-cache.
    """
    digest = hashlib.md5(f"{peptidoform}|{position}".encode()).digest()
    return int.from_bytes(digest[:4], "big", signed=False) & 0x7FFFFFFF


_PICKER_REF_SCHEMA = {
    "id": pl.Int32, "mzml": pl.String, "scan": pl.Int32,
    "rt": pl.Float32, "charge": pl.Int32, "score": pl.Float32,
}
_PICKER_QUERY_SCHEMA = {
    "psm_id": pl.Int32, "mzml": pl.String, "scan": pl.Int32,
    "rt": pl.Float32, "charge": pl.Int32, "score": pl.Float32,
}


def _build_picker_frames(
    peptidoform: str, annotated_queries: dict, reference_store: dict, psm_id_to_peptidoform: dict,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Build (reference, query) LazyFrames for the picker — refs use hashed ids, queries use real psm_id.
    """
    def picker_columns(psm) -> dict:
        # Each spectrum is a carrier PSM; read the display fields off it.
        # score is None for no search score (Koina) or NaN
        raw = psm.score
        score = float(raw) if raw is not None and raw == raw else None
        return {
            "mzml": source_label(psm),
            "scan": extract_id_number(psm.spectrum_id),
            "rt": float(psm.retention_time) if psm.retention_time is not None else None,
            "charge": int(psm.get_precursor_charge() or 0),
            "score": score,
        }

    # Reference rows: hashed per (peptidoform, position) — see _ref_row_id docstring
    ref_rows = [
        {"id": _ref_row_id(peptidoform, i), **picker_columns(psm)}
        for i, (psm, _) in enumerate(reference_store.get(peptidoform, []))
    ]
    # Query rows: filter all annotated queries down to the chosen peptidoform
    query_rows = [
        {"psm_id": psm_id, **picker_columns(psm)}
        for psm_id, (psm, _) in annotated_queries.items()
        if psm_id_to_peptidoform.get(psm_id) == peptidoform
    ]

    # Materialise: typed-empty LazyFrame when no rows so the Table still renders
    ref_frame = (
        pl.DataFrame(ref_rows, schema=_PICKER_REF_SCHEMA).lazy() if ref_rows
        else pl.DataFrame(schema=_PICKER_REF_SCHEMA).lazy()
    )
    query_frame = (
        pl.DataFrame(query_rows, schema=_PICKER_QUERY_SCHEMA).lazy() if query_rows
        else pl.DataFrame(schema=_PICKER_QUERY_SCHEMA).lazy()
    )
    return ref_frame, query_frame


def _picker_column_definitions() -> list[dict]:
    """Return Tabulator column defs for the picker — same money formatter idiom as explorer.py.
    """
    return [
        {"field": "mzml", "title": "mzML", "sorter": "string"},
        {"field": "scan", "title": "Scan", "sorter": "number", "hozAlign": "right"},
        {"field": "rt", "title": "RT", "sorter": "number",
         "formatter": "money",
         "formatterParams": {"precision": 2, "symbol": ""},
         "hozAlign": "right"},
        {"field": "charge", "title": "Z", "sorter": "number", "hozAlign": "right"},
        {"field": "score", "title": "Score", "sorter": "number",
         "formatter": "money",
         "formatterParams": {"precision": 4, "symbol": ""},
         "hozAlign": "right"},
    ]


def _resolve_picker_cache_path(params: dict) -> str:
    """Return the picker's component-cache directory (same dir as Explorer, distinct cache_ids).
    """
    export_dir = get_param(params, "export_directory", "")
    cache_path = export_dir or str(Path(st.session_state.workspace) / "exports")
    os.makedirs(cache_path, exist_ok=True)
    return cache_path


def _get_picker_state_manager() -> StateManager:
    """Lazily build the picker's StateManager (distinct session_key from the Explorer).
    """
    if "mspecviz_mirror_picker_sm" not in st.session_state:
        st.session_state["mspecviz_mirror_picker_sm"] = StateManager(
            session_key="mspecviz_mirror_picker_state",
        )
    return st.session_state["mspecviz_mirror_picker_sm"]


def _clear_picker_selections(state_manager: StateManager) -> None:
    """Reset both StateManager selections AND Tabulator's component-local trackers.
    regenerate_cache=True only rewrites the parquet cache; without this the old row stays highlighted.
    """
    state_manager.clear_selection("selected_ref_row")
    state_manager.clear_selection("selected_query_psm")
    for tracker_key in ("_svc_table_last_selection", "_svc_table_last_sort_filter"):
        cache_dict = st.session_state.get(tracker_key, {})
        cache_dict.pop("mirror_ref_picker", None)
        cache_dict.pop("mirror_query_picker", None)


def _render_metrics_expander(metrics: dict) -> None:
    """Show the 7 similarity metrics + ion counts for one query/reference pair in a collapsed expander.
    """
    with st.expander("Comparison Metrics", expanded=False):
        metric_data = {
            METRIC_SHORT_TO_DISPLAY[key]: (
                f"{value:.4f}" if not (isinstance(value, float) and value != value) else "N/A"
            )
            for key, value in metrics.items() if key in METRIC_SHORT_TO_DISPLAY
        }
        metric_data["Shared Ions"] = str(metrics.get("shared", 0))
        metric_data["Annotated (Query)"] = str(metrics.get("query_annotation_count", 0))
        metric_data["Annotated (Ref)"] = str(metrics.get("reference_annotation_count", 0))
        st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)


def _render_mirror_pair(query_psm, query_frame,
                        reference_psm, reference_frame) -> None:
    """Compute metrics, draw the Plotly mirror plot, and show the metrics expander below it.
    """
    metrics = compute_all_metrics(query_frame, reference_frame)
    query_label, ref_label = format_spectrum_labels(query_psm, reference_psm)
    title = f"{query_psm.peptidoform.modified_sequence} +{query_psm.get_precursor_charge() or 0}"
    fig = create_spectrum_plot(
        query_frame, reference_frame,
        backend="ms_plotly",
        title=title, query_label=query_label, reference_label=ref_label,
        metrics_text=format_metrics_text(metrics),
    )
    st.plotly_chart(fig, width="stretch")
    _render_metrics_expander(metrics)


def _render_picker_tables(
    ref_frame: pl.LazyFrame, query_frame: pl.LazyFrame,
    cache_path: str, state_manager: StateManager,
) -> None:
    """Render the two side-by-side openms_insight Tables wired to *state_manager*.
    """
    columns = _picker_column_definitions()
    left_column, right_column = st.columns(2)
    # default_row=-1 disables Tabulator's row-0 auto-selection on first render;
    # without it Table eagerly fires interactivity for row 0 and the "both
    # picks present" gate below is bypassed. Experimental on the left matches
    # the mirror-plot convention (query = upper trace).
    with left_column:
        st.subheader("Experimental")
        st.caption("Upper spectrum in the mirror plot below.")
        Table(
            cache_id="mirror_query_picker", data=query_frame, cache_path=cache_path,
            column_definitions=columns, index_field="psm_id", regenerate_cache=True,
            default_row=-1,
            interactivity={"selected_query_psm": "psm_id"},
        )(key="mirror_query_picker", state_manager=state_manager)
    with right_column:
        st.subheader("Reference")
        st.caption("Lower spectrum in the mirror plot below.")
        Table(
            cache_id="mirror_ref_picker", data=ref_frame, cache_path=cache_path,
            column_definitions=columns, index_field="id", regenerate_cache=True,
            default_row=-1,
            interactivity={"selected_ref_row": "id"},
        )(key="mirror_ref_picker", state_manager=state_manager)


@st.fragment
def mirror_plots_fragment(params: dict):
    """Render the Mirror Plots sub-view — peptide selectbox -> two pickers -> mirror plot.
    """
    # Pull session state: annotated queries, references, and peptidoform PSMs from Load & Annotate
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    reference_store = st.session_state.get(get_state_key("reference_store"), {})
    peptidoform_psms = st.session_state.get(get_state_key("peptidoform_psms"), {})

    if not annotated_queries:
        st.info("No annotated spectra available.")
        return

    # Eligibility: peptidoforms needing both a query PSM and a reference to mirror.
    # Reshape the session-state dict into the canonical [(psm_id, query_psm, query_frame)]
    # list the core mapping helper expects.
    annotated_query_list = [
        (psm_id, query_psm, query_frame)
        for psm_id, (query_psm, query_frame) in annotated_queries.items()
    ]
    psm_id_to_peptidoform = build_psm_to_peptidoform_mapping(peptidoform_psms, annotated_query_list)
    eligible = _eligible_peptidoforms(peptidoform_psms, reference_store, psm_id_to_peptidoform)

    # Hidden caption: surface how many peptidoforms were dropped for missing a reference pair
    peptidoforms_with_query = set(psm_id_to_peptidoform.values())
    peptidoforms_with_query_count = sum(1 for peptidoform in peptidoform_psms if peptidoform in peptidoforms_with_query)
    hidden_peptidoform_count = peptidoforms_with_query_count - len(eligible)
    if hidden_peptidoform_count:
        st.caption(
            f"{hidden_peptidoform_count} peptide(s) hidden — no reference / experimental pair available."
        )
    if not eligible:
        query_count = len(annotated_queries)
        reference_count = sum(len(refs) for refs in reference_store.values())
        st.info(
            f"No peptides have both a reference and an experimental spectrum "
            f"({query_count} queries / {reference_count} references loaded)."
        )
        return

    # Peptide selector: snap a stale selection back to a still-eligible peptidoform
    selectbox_key = f"{STATE_PREFIX}_mirror_chosen_peptidoform"
    if st.session_state.get(selectbox_key) not in eligible:
        st.session_state[selectbox_key] = eligible[0]
    chosen = st.selectbox(
        "Peptide", options=eligible,
        format_func=lambda peptidoform: peptidoform,
        key=selectbox_key,
    )

    # Peptidoform switch: clear Tabulator + StateManager state so stale picks don't carry over
    state_manager = _get_picker_state_manager()
    last_peptidoform = f"{STATE_PREFIX}_mirror_last_peptidoform"
    if st.session_state.get(last_peptidoform) != chosen:
        _clear_picker_selections(state_manager)
        st.session_state[last_peptidoform] = chosen

    # Picker tables: two side-by-side openms_insight Tables for the chosen peptidoform
    reference_picker_frame, query_picker_frame = _build_picker_frames(
        chosen, annotated_queries, reference_store, psm_id_to_peptidoform,
    )
    _render_picker_tables(
        reference_picker_frame, query_picker_frame, _resolve_picker_cache_path(params), state_manager,
    )

    # Read selections from the StateManager. openms_insight's bridge auto-fills
    # both keys with row 0 on first render (bridge.py:694-697, unconditional —
    # ignores default_row), so the plot renders for the first ref / first
    # query of the chosen peptide; clicking either table swaps just that
    # side. None-checks below are defensive for AppTest (no Vue running).
    state = state_manager.get_state_for_vue()
    ref_id = state.get("selected_ref_row")
    selected_query_psm_id = state.get("selected_query_psm")
    if ref_id is None or selected_query_psm_id is None:
        st.info("Pick one row in each table to render the mirror plot.")
        return

    # Resolve picks: map both row ids back to the in-memory spectra
    if selected_query_psm_id not in annotated_queries:
        return
    refs = reference_store.get(chosen, [])
    # Rebuild the hashed-id reverse index for this peptidoform so a stale ref_id
    # from a previous peptide (after auto-fill races) resolves to None and we
    # bail safely. See _ref_row_id for the hashing rationale.
    ref_index = {_ref_row_id(chosen, i): i for i in range(len(refs))}
    if ref_id not in ref_index:
        return

    # Render: compute metrics and draw the mirror plot for the chosen pair
    query_psm, query_frame = annotated_queries[selected_query_psm_id]
    reference_psm, reference_frame = refs[ref_index[ref_id]]
    _render_mirror_pair(query_psm, query_frame, reference_psm, reference_frame)
