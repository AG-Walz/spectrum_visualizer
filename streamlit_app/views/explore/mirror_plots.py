"""
Interactive Plotly mirror-plot fragment.

Two-stage picker: a peptide selectbox over variants that have BOTH a reference
and an experimental spectrum loaded, then two side-by-side openms_insight
Tables (left = references, right = experimentals). The mirror plot renders
only when both sides have a single-row selection. Linked-selection style
mirrors the Explorer fragment (toppview-lite reference).
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
from mspecviz.annotation import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    format_metrics_text,
)
from mspecviz.common.defaults import DEFAULTS
from mspecviz.pipeline import (
    build_psm_to_variant_mapping,
    format_spectrum_labels,
)
from mspecviz.plotting import create_plotly_mirror_plot
from mspecviz.spectra.mzml_native_ids import _extract_id_number


def _eligible_variants(variant_groups: dict, reference_store: dict, psm_to_vk: dict) -> list:
    """Variants with at least one query AND one reference loaded."""
    has_query = set(psm_to_vk.values())
    return [vk for vk in variant_groups
            if vk in has_query and reference_store.get(vk)]


def _picker_row(spec, score) -> dict:
    """Translate one ``LoadedSpectrum`` into the columns the Tables display."""
    return {
        "mzml": spec.source_label,
        "scan": _extract_id_number(spec.spectrum_id),
        "rt": float(spec.rt) if spec.rt is not None else None,
        "charge": int(spec.precursor_charge or 0),
        "score": score,
    }


def _safe_score(spec) -> float | None:
    """Pull ``match.score`` if present; return ``None`` for NaN or no-match (Koina)."""
    if spec.match is None:
        return None
    score = spec.match.score
    return float(score) if score == score else None  # NaN check


def _ref_row_id(variant_key: str, position: int) -> int:
    """Deterministic, peptide-scoped row id for the reference picker.

    Positional ids ([0, 1, 2, …]) collide across peptides, which makes the
    Vue render cache key for the reference Table identical between peptide
    switches (the first auto-filled selection is always 0) and the table
    therefore appears stuck on stale data. Hashing in the *variant_key*
    makes the auto-filled ``selected_ref_row`` differ per peptide, busting
    the cache and forcing the bridge to push fresh rows. The query Table
    doesn't need this — its index is the real ``psm_id`` which is already
    unique across peptides.
    """
    digest = hashlib.md5(f"{variant_key}|{position}".encode()).digest()
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
    variant_key: str, annotated_queries: dict, reference_store: dict, psm_to_vk: dict,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Two LazyFrames — references (id = positional) and queries (psm_id = real key)."""
    ref_rows = [
        {"id": _ref_row_id(variant_key, i), **_picker_row(spec, _safe_score(spec))}
        for i, (spec, _) in enumerate(reference_store.get(variant_key, []))
    ]
    query_rows = [
        {"psm_id": psm_id, **_picker_row(spec, _safe_score(spec))}
        for psm_id, (spec, _) in annotated_queries.items()
        if psm_to_vk.get(psm_id) == variant_key
    ]

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
    """Tabulator column defs — same ``money`` formatter idiom as ``explorer.py``."""
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
    """Mirror of ``explorer._resolve_cache_path`` — same dir, distinct cache_ids."""
    export_dir = get_param(params, "export_directory", "")
    cache_path = export_dir or str(Path(st.session_state.workspace) / "exports")
    os.makedirs(cache_path, exist_ok=True)
    return cache_path


def _get_picker_state_manager() -> StateManager:
    """Lazy StateManager for the picker; distinct session_key from the Explorer."""
    if "mspecviz_mirror_picker_sm" not in st.session_state:
        st.session_state["mspecviz_mirror_picker_sm"] = StateManager(
            session_key="mspecviz_mirror_picker_state",
        )
    return st.session_state["mspecviz_mirror_picker_sm"]


def _clear_picker_selections(sm: StateManager) -> None:
    """Reset both StateManager keys AND Tabulator's component-local trackers.

    ``regenerate_cache=True`` only rewrites the parquet cache — Tabulator's
    ``_LAST_SELECTION_KEY`` / ``_LAST_SORT_FILTER_KEY`` survive it, so we pop
    them explicitly here. Without this the previously-clicked row would stay
    visually highlighted after a peptide switch.
    """
    sm.clear_selection("selected_ref_row")
    sm.clear_selection("selected_query_psm")
    for tracker in ("_svc_table_last_selection", "_svc_table_last_sort_filter"):
        d = st.session_state.get(tracker, {})
        d.pop("mirror_ref_picker", None)
        d.pop("mirror_query_picker", None)


def _render_metrics_expander(metrics: dict) -> None:
    """Render the per-comparison metrics table inside an expander."""
    with st.expander("Comparison Metrics", expanded=False):
        metric_data = {
            METRIC_SHORT_TO_DISPLAY[k]: f"{v:.4f}" if not (isinstance(v, float) and v != v) else "N/A"
            for k, v in metrics.items() if k in METRIC_SHORT_TO_DISPLAY
        }
        metric_data["Shared Ions"] = str(metrics.get("shared", 0))
        metric_data["Annotated (Query)"] = str(metrics.get("ann_q", 0))
        metric_data["Annotated (Ref)"] = str(metrics.get("ann_r", 0))
        st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)


def _render_mirror_pair(query_spectrum, query_frame,
                        reference_spectrum, reference_frame, normalize: bool) -> None:
    """Render one mirror plot + metrics expander for the chosen query/reference pair."""
    metrics = compute_all_metrics(query_frame, reference_frame)
    query_label, ref_label = format_spectrum_labels(query_spectrum, reference_spectrum)
    title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"
    fig = create_plotly_mirror_plot(
        upper_spectrum=query_frame, lower_spectrum=reference_frame,
        title=title, upper_label=query_label, lower_label=ref_label,
        metrics_text=format_metrics_text(metrics), normalize=normalize,
    )
    st.plotly_chart(fig, width="stretch")
    _render_metrics_expander(metrics)


def _render_picker_tables(
    ref_frame: pl.LazyFrame, query_frame: pl.LazyFrame,
    cache_path: str, sm: StateManager,
) -> None:
    """Render the two side-by-side openms_insight Tables wired to *sm*."""
    cols = _picker_column_definitions()
    col_left, col_right = st.columns(2)
    # default_row=-1 disables Tabulator's auto-selection on first render; without
    # it the gate "render plot only when both selections exist" is bypassed
    # immediately because Table eagerly fires interactivity for row 0.
    with col_left:
        st.subheader("Reference")
        Table(
            cache_id="mirror_ref_picker", data=ref_frame, cache_path=cache_path,
            column_definitions=cols, index_field="id", regenerate_cache=True,
            default_row=-1,
            interactivity={"selected_ref_row": "id"},
        )(key="mirror_ref_picker", state_manager=sm)
    with col_right:
        st.subheader("Experimental")
        Table(
            cache_id="mirror_query_picker", data=query_frame, cache_path=cache_path,
            column_definitions=cols, index_field="psm_id", regenerate_cache=True,
            default_row=-1,
            interactivity={"selected_query_psm": "psm_id"},
        )(key="mirror_query_picker", state_manager=sm)


@st.fragment
def mirror_plots_fragment(params: dict):
    """Two-stage picker → side-by-side selection → mirror plot."""
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    reference_store = st.session_state.get(get_state_key("reference_store"), {})
    variant_groups = st.session_state.get(get_state_key("variant_groups"), {})

    if not annotated_queries:
        st.info("No annotated spectra available.")
        return

    psm_to_vk = build_psm_to_variant_mapping(variant_groups, annotated_queries)
    eligible = _eligible_variants(variant_groups, reference_store, psm_to_vk)

    has_query = set(psm_to_vk.values())
    queryful = sum(1 for vk in variant_groups if vk in has_query)
    hidden = queryful - len(eligible)
    if hidden:
        st.caption(
            f"{hidden} peptide(s) hidden — no reference / experimental pair available."
        )
    if not eligible:
        n_q = len(annotated_queries)
        n_r = sum(len(v) for v in reference_store.values())
        st.info(
            f"No peptides have both a reference and an experimental spectrum "
            f"({n_q} queries / {n_r} references loaded)."
        )
        return

    sel_key = f"{STATE_PREFIX}_mirror_chosen_variant"
    if st.session_state.get(sel_key) not in eligible:
        st.session_state[sel_key] = eligible[0]
    chosen = st.selectbox(
        "Peptide", options=eligible,
        format_func=lambda vk: variant_groups[vk].display_label,
        key=sel_key,
    )

    sm = _get_picker_state_manager()
    last_key = f"{STATE_PREFIX}_mirror_last_variant"
    if st.session_state.get(last_key) != chosen:
        _clear_picker_selections(sm)
        st.session_state[last_key] = chosen

    ref_frame, query_frame = _build_picker_frames(
        chosen, annotated_queries, reference_store, psm_to_vk,
    )
    _render_picker_tables(ref_frame, query_frame, _resolve_picker_cache_path(params), sm)

    # openms_insight's bridge auto-fills both selection keys with the first
    # row's value on first render (bridge.py:694-697 — runs unconditionally
    # regardless of default_row). So the plot renders immediately for the
    # first ref / first query of the chosen peptide; clicking either table
    # swaps just that side. The None-checks below are defensive — Vue isn't
    # running under AppTest, so the auto-fill may not fire there.
    state = sm.get_state_for_vue()
    ref_id = state.get("selected_ref_row")
    query_psm = state.get("selected_query_psm")
    if ref_id is None or query_psm is None:
        st.info("Pick one row in each table to render the mirror plot.")
        return

    if query_psm not in annotated_queries:
        return
    refs = reference_store.get(chosen, [])
    # ref ids are hashed per (variant_key, position) — see _ref_row_id docstring.
    # Rebuild the reverse index here so a stale ref_id from a previous peptide
    # (after auto-fill races) resolves to None and we bail safely.
    ref_index = {_ref_row_id(chosen, i): i for i in range(len(refs))}
    if ref_id not in ref_index:
        return

    query_spectrum, query_anno = annotated_queries[query_psm]
    reference_spectrum, reference_anno = refs[ref_index[ref_id]]
    normalize = get_param(params, "intensity_mode", DEFAULTS.intensity_mode) == "Relative"
    _render_mirror_pair(
        query_spectrum, query_anno, reference_spectrum, reference_anno, normalize,
    )
