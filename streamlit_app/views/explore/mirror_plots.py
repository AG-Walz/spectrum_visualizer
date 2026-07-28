"""Mirror Plots sub-view of the Explore tab — Plotly mirror plot driven by a two-stage picker.
Rendered by explore_page.py: peptide selectbox + two Tables (query left, reference right).
"""
from __future__ import annotations

import hashlib
import io
import os
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st
from openms_insight import StateManager, Table
from psm_utils import Peptidoform

from streamlit_app.utils.state import STATE_PREFIX, get_param, get_state_key
from core.comparison_metrics import (
    METRIC_KEY_TO_COLUMN_NAME,
    compute_all_metrics,
    format_metrics_text,
)
from core.defaults import DEFAULTS
from core.export import sanitize_filename
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


# Similarity metrics the best-spectrum auto-pick can rank by — these labels are the
# comparison_frame column names, all normalised so higher = more similar.
_SIMILARITY_METRICS = list(METRIC_KEY_TO_COLUMN_NAME.values())


def _modification_name(modification) -> str:
    """Best-effort display name for a psm_utils modification object."""
    return getattr(modification, "name", None) or getattr(modification, "value", None) or str(modification)


def _modification_descriptor(peptidoform: Peptidoform) -> str:
    """Compact human label for a peptidoform's modifications ("N-term(Acetyl), M3(Oxidation)").
    Returns "(unmodified)" when the peptidoform carries no modification.
    """
    if not peptidoform.is_modified:
        return "(unmodified)"
    parts = [f"N-term({_modification_name(m)})" for m in (peptidoform.properties.get("n_term") or [])]
    for position, (amino_acid, modifications) in enumerate(peptidoform.parsed_sequence, start=1):
        parts.extend(f"{amino_acid}{position}({_modification_name(m)})" for m in (modifications or []))
    parts.extend(f"C-term({_modification_name(m)})" for m in (peptidoform.properties.get("c_term") or []))
    return ", ".join(parts) if parts else "(unmodified)"


def _build_peptidoform_index(eligible: list) -> tuple[dict, dict]:
    """Split eligible peptidoform keys into a sequence -> [keys] map and a key -> descriptor map.
    The peptidoform key is a modified-sequence string (no charge); base sequence and modification
    label are parsed from it. Keys under each sequence are sorted by their modification descriptor;
    sequences are sorted alphabetically by the caller (#106).
    """
    sequence_to_keys: dict[str, list[str]] = {}
    descriptor_of: dict[str, str] = {}
    for key in eligible:
        peptidoform = Peptidoform(key)
        sequence_to_keys.setdefault(peptidoform.sequence.upper(), []).append(key)
        descriptor_of[key] = _modification_descriptor(peptidoform)
    for keys in sequence_to_keys.values():
        keys.sort(key=lambda key: descriptor_of[key])
    return sequence_to_keys, descriptor_of


def _query_charges(peptidoform: str, annotated_queries: dict, psm_id_to_peptidoform: dict) -> list[int]:
    """Return the sorted distinct precursor charges among a peptidoform's query spectra."""
    return sorted({
        int(query_psm.get_precursor_charge() or 0)
        for psm_id, (query_psm, _frame) in annotated_queries.items()
        if psm_id_to_peptidoform.get(psm_id) == peptidoform
    })


def _best_pair_by_metric(
    comparison_frame, peptidoform_psm_ids: set, reference_store: dict,
    peptidoform: str, charges: list[int], metric_column: str,
) -> tuple[float, int, int] | None:
    """Find the (query, reference) pair with the highest similarity for one peptidoform.
    Reads the already-computed metric straight from comparison_frame (no recompute) and resolves the
    winning row's ref_spectrum_id back to a reference_store index. Returns
    (metric_value, query_psm_id, reference_index), or None when nothing scores. metric_column is any
    comparison_frame similarity column (e.g. "R² match", "Cosine"); all are normalised higher = better.
    """
    if comparison_frame is None or comparison_frame.empty or metric_column not in comparison_frame.columns:
        return None
    rows = comparison_frame[comparison_frame["psm_id"].isin(peptidoform_psm_ids)]
    if charges:
        rows = rows[rows["Charge"].isin(charges)]
    rows = rows.dropna(subset=[metric_column])
    if rows.empty:
        return None
    best = rows.loc[rows[metric_column].idxmax()]
    # Resolve on (source label, spectrum_id): spectrum_id alone is not unique across mzML files,
    # so a bare-id map would collapse colliding native IDs to the wrong reference.
    reference_index = {(source_label(reference_psm), reference_psm.spectrum_id): index
                       for index, (reference_psm, _frame) in enumerate(reference_store.get(peptidoform, []))
                       }.get((best["Reference"], best["ref_spectrum_id"]))
    if reference_index is None:
        return None
    return float(best[metric_column]), int(best["psm_id"]), reference_index


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
    charges: list[int] | None = None,
) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Build (reference, query) LazyFrames for the picker — refs use hashed ids, queries use real psm_id.
    charges (when non-empty) restricts the query rows to those precursor charges; references are
    always listed in full so a match at any charge stays reachable.
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
    # Query rows: filter all annotated queries down to the chosen peptidoform (and charge, if set)
    query_rows = [
        {"psm_id": psm_id, **picker_columns(psm)}
        for psm_id, (psm, _) in annotated_queries.items()
        if psm_id_to_peptidoform.get(psm_id) == peptidoform
        and (not charges or int(psm.get_precursor_charge() or 0) in charges)
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
            METRIC_KEY_TO_COLUMN_NAME[key]: (
                f"{value:.4f}" if not (isinstance(value, float) and value != value) else "N/A"
            )
            for key, value in metrics.items() if key in METRIC_KEY_TO_COLUMN_NAME
        }
        metric_data["Shared Ions"] = str(metrics.get("shared", 0))
        metric_data["Annotated (Query)"] = str(metrics.get("query_annotation_count", 0))
        metric_data["Annotated (Ref)"] = str(metrics.get("reference_annotation_count", 0))
        st.dataframe(pd.DataFrame([metric_data]), width="stretch", hide_index=True)


def _static_mirror_bytes(query_frame: pd.DataFrame, reference_frame: pd.DataFrame, *,
                         title: str, query_label: str, reference_label: str,
                         metrics_text: str, max_peak_labels: int,
                         export_resolution: int) -> tuple[bytes, bytes]:
    """Render the static (matplotlib) twin of the on-screen plot and return its (PNG, SVG) bytes.
    create_spectrum_plot builds a bare Figure (no pyplot registry), so it needs no explicit close.
    """
    fig = create_spectrum_plot(
        query_frame, reference_frame,
        backend="ms_matplotlib",
        title=title, query_label=query_label, reference_label=reference_label,
        metrics_text=metrics_text, max_peak_labels=max_peak_labels,
    )
    rendered = []
    for output_format in ("png", "svg"):
        buffer = io.BytesIO()
        # DPI only meaningfully applies to the raster PNG; SVG is vector.
        save_kwargs = {"dpi": export_resolution} if output_format == "png" else {}
        fig.savefig(buffer, format=output_format, bbox_inches="tight", facecolor="white",
                    **save_kwargs)
        rendered.append(buffer.getvalue())
    return rendered[0], rendered[1]


# Cached view of the renderer above, so clicking PNG and then SVG for the same pair renders
# once. The undecorated function stays importable for tests.
_cached_static_mirror_bytes = st.cache_data(show_spinner=False, max_entries=8)(_static_mirror_bytes)


def _download_basename(query_psm, reference_psm) -> str:
    """Filesystem-safe stem for a downloaded mirror plot, matching the Export tab's naming.
    """
    charge = query_psm.get_precursor_charge() or 0
    return sanitize_filename(
        f"{query_psm.peptidoform.modified_sequence}_z{charge}"
        f"_{source_label(query_psm)}_{query_psm.spectrum_id}"
        f"_vs_{source_label(reference_psm)}_{reference_psm.spectrum_id}"
    )


def _render_download_row(fig, query_frame, reference_frame, *, basename: str,
                         title: str, query_label: str, reference_label: str,
                         metrics_text: str, max_peak_labels: int, export_resolution: int) -> None:
    """Offer the plot on screen as PNG / SVG / standalone HTML downloads (#107).
    PNG and SVG come from the matplotlib backend so a saved plot matches what the Export tab
    writes; the HTML is the interactive Plotly figure already on screen, pulling plotly.js from
    the CDN exactly like the exported HTML plots do. Each payload is passed as a callable so
    nothing is rendered until the user actually clicks — same lazy idiom as the export ZIP.
    """
    def static_bytes(index: int):
        return lambda: _cached_static_mirror_bytes(
            query_frame, reference_frame, title=title, query_label=query_label,
            reference_label=reference_label, metrics_text=metrics_text,
            max_peak_labels=max_peak_labels, export_resolution=export_resolution,
        )[index]

    st.caption("Download this plot")
    png_col, svg_col, html_col = st.columns(3)
    with png_col:
        st.download_button("PNG", data=static_bytes(0), file_name=f"{basename}.png",
                           mime="image/png", icon=":material/download:",
                           key=f"{STATE_PREFIX}_mirror_download_png")
    with svg_col:
        st.download_button("SVG", data=static_bytes(1), file_name=f"{basename}.svg",
                           mime="image/svg+xml", icon=":material/download:",
                           key=f"{STATE_PREFIX}_mirror_download_svg")
    with html_col:
        st.download_button("Interactive HTML",
                           data=lambda: fig.to_html(include_plotlyjs="cdn"),
                           file_name=f"{basename}.html", mime="text/html",
                           icon=":material/download:",
                           key=f"{STATE_PREFIX}_mirror_download_html")


def _render_mirror_pair(query_psm, query_frame,
                        reference_psm, reference_frame,
                        max_peak_labels: int = 0,
                        export_resolution: int = DEFAULTS.export_resolution) -> None:
    """Compute metrics, draw the Plotly mirror plot, then the download row and metrics expander.
    """
    metrics = compute_all_metrics(query_frame, reference_frame)
    query_label, ref_label = format_spectrum_labels(query_psm, reference_psm)
    title = f"{query_psm.peptidoform.modified_sequence} +{query_psm.get_precursor_charge() or 0}"
    metrics_text = format_metrics_text(metrics)
    fig = create_spectrum_plot(
        query_frame, reference_frame,
        backend="ms_plotly",
        title=title, query_label=query_label, reference_label=ref_label,
        metrics_text=metrics_text,
        max_peak_labels=max_peak_labels,
    )
    st.plotly_chart(fig, width="stretch")
    _render_download_row(
        fig, query_frame, reference_frame,
        basename=_download_basename(query_psm, reference_psm),
        title=title, query_label=query_label, reference_label=ref_label,
        metrics_text=metrics_text, max_peak_labels=max_peak_labels,
        export_resolution=export_resolution,
    )
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
    """Render the Mirror Plots sub-view — sequence/modification/charge filters -> mirror plot.
    Either auto-picks the highest-R² query/reference pair or drives two manual pickers (#62/#106).
    """
    # Pull session state: annotated queries, references, and peptidoform PSMs from Load & Annotate
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    reference_store = st.session_state.get(get_state_key("reference_store"), {})
    peptidoform_psms = st.session_state.get(get_state_key("peptidoform_psms"), {})
    comparison_frame = st.session_state.get(get_state_key("comparison_frame"))

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

    sequence_to_keys, descriptor_of = _build_peptidoform_index(eligible)
    sequences = sorted(sequence_to_keys)

    # Row 1 — WHICH SPECTRA: pick the peptide (sequence + modification, #62/#106) and charge.
    sequence_col, modification_col, charge_col = st.columns([3, 3, 2])
    with sequence_col:
        sequence_key = f"{STATE_PREFIX}_mirror_sequence"
        if st.session_state.get(sequence_key) not in sequences:
            st.session_state[sequence_key] = sequences[0]
        chosen_sequence = st.selectbox("Peptide sequence", options=sequences, key=sequence_key)
    with modification_col:
        modifications = sequence_to_keys[chosen_sequence]
        modification_key = f"{STATE_PREFIX}_mirror_modification"
        if st.session_state.get(modification_key) not in modifications:
            st.session_state[modification_key] = modifications[0]
        chosen = st.selectbox("Modification", options=modifications,
                              format_func=descriptor_of.get, key=modification_key)

    # Peptidoform switch: clear stale picks and re-default the charge filter to all.
    state_manager = _get_picker_state_manager()
    last_peptidoform = f"{STATE_PREFIX}_mirror_last_peptidoform"
    if st.session_state.get(last_peptidoform) != chosen:
        _clear_picker_selections(state_manager)
        st.session_state[last_peptidoform] = chosen
        st.session_state.pop(f"{STATE_PREFIX}_mirror_charges", None)

    available_charges = _query_charges(chosen, annotated_queries, psm_id_to_peptidoform)
    charges_key = f"{STATE_PREFIX}_mirror_charges"
    st.session_state[charges_key] = [c for c in st.session_state.get(charges_key, available_charges)
                                     if c in available_charges]  # drop charges stale after a reload
    with charge_col:
        selected_charges = st.multiselect(
            "Charge", options=available_charges, format_func=lambda charge: f"+{charge}",
            key=charges_key, help="Filter the experimental spectra by precursor charge (empty = all).",
        )

    # Row 2 — HOW TO PLOT: label density, then the best-spectrum auto-pick and the metric it ranks by.
    # bottom-aligned so the checkbox lines up with the labelled widgets beside it.
    labels_col, auto_col, measure_col = st.columns([2, 3, 2], vertical_alignment="bottom")
    with labels_col:
        # Namespaced key keeps this an independent live override of the persisted export setting.
        max_peak_labels = st.number_input(
            "Max peak labels", min_value=0, max_value=200, step=5,
            value=int(get_param(params, "max_peak_labels", DEFAULTS.max_peak_labels)),
            help="Label only the N most intense annotated peaks per spectrum (0 = all).",
            key=f"{STATE_PREFIX}_mirror_max_peak_labels",
        )
    with auto_col:
        auto_best = st.checkbox("Show only the best spectrum by similarity metric",
                                key=f"{STATE_PREFIX}_mirror_auto_best")
    with measure_col:
        best_metric = st.selectbox("Similarity metric", options=_SIMILARITY_METRICS,
                                   disabled=not auto_best, key=f"{STATE_PREFIX}_mirror_best_metric")

    # Plots downloaded from this view reuse the workspace's export resolution, so a PNG saved
    # here matches one written by the Export tab (#107).
    export_resolution = int(get_param(params, "export_resolution", DEFAULTS.export_resolution))

    # Auto-pick mode: skip the tables and render the highest-scoring pair directly.
    if auto_best:
        peptidoform_psm_ids = {psm_id for psm_id, peptidoform in psm_id_to_peptidoform.items()
                               if peptidoform == chosen}
        best = _best_pair_by_metric(comparison_frame, peptidoform_psm_ids, reference_store,
                                    chosen, selected_charges, best_metric)
        if best is None:
            st.info(f"No query/reference pair has a computable {best_metric} for this selection.")
            return
        metric_value, best_psm_id, best_reference_index = best
        st.caption(f"Best pair by {best_metric}: {metric_value:.4f}")
        query_psm, query_frame = annotated_queries[best_psm_id]
        reference_psm, reference_frame = reference_store[chosen][best_reference_index]
        _render_mirror_pair(query_psm, query_frame, reference_psm, reference_frame,
                            max_peak_labels, export_resolution)
        return

    # Manual mode: two side-by-side pickers (query rows filtered by charge).
    reference_picker_frame, query_picker_frame = _build_picker_frames(
        chosen, annotated_queries, reference_store, psm_id_to_peptidoform, selected_charges,
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
    _render_mirror_pair(query_psm, query_frame, reference_psm, reference_frame,
                        max_peak_labels, export_resolution)
