"""
Shared processing logic for the MSInsight pipeline.

Used by both ``cli.py`` (terminal mode) and ``content/peptide_mirror.py``
(Streamlit mode).  Keeps all duplicated spectrum-processing helpers in a
single authoritative location.
"""
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

logger = logging.getLogger("msinsight.pipeline")


# ============================================================================
# Label helpers
# ============================================================================

def get_source_label(spec) -> str:
    """Extract filename stem from spectrum, or 'Koina' if no source file."""
    return Path(spec.source_file).stem if spec.source_file else "Koina"


def format_spectrum_labels(query_spec, ref_spec=None):
    """Format query/reference spectrum labels with RT.

    Returns ``(query_label, ref_label)`` tuple.  *ref_label* is ``""``
    when *ref_spec* is ``None``.
    """
    query_label = f"Query: {get_source_label(query_spec)} (RT: {query_spec.rt:.1f})"
    if ref_spec is None:
        return query_label, ""
    ref_label = f"Ref: {get_source_label(ref_spec)} (RT: {ref_spec.rt:.1f})"
    return query_label, ref_label


# ============================================================================
# PSM-to-variant mapping
# ============================================================================

def build_psm_to_variant_mapping(variant_groups, annotated_queries):
    """Build mapping from PSM ID to variant key.

    *annotated_queries* can be:

    * ``dict``  -- ``{psm_id: (spec, df)}``   (Streamlit format)
    * ``list``  -- ``[(psm_id, spec, df)]``    (CLI format)
    """
    from src.spectrum_loader import spectrum_lookup_key

    psm_id_to_vk: Dict[str, str] = {}

    # Normalise to iterable of (pid, spec)
    if isinstance(annotated_queries, dict):
        items = [(pid, spec) for pid, (spec, _) in annotated_queries.items()]
    else:
        items = [(pid, spec) for pid, spec, _ in annotated_queries]

    for vk, vg in variant_groups.items():
        for match in vg.matches:
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for pid, spec in items:
                if spectrum_lookup_key(spec) == key:
                    psm_id_to_vk[pid] = vk
                    break
    return psm_id_to_vk


# ============================================================================
# Koina prediction
# ============================================================================

def get_koina_prediction(peptidoform, charge, annotation_kwargs,
                         koina_model="ms2pip_timsTOF2024",
                         koina_ce=35, koina_instrument="TIMSTOF"):
    """Get Koina AI prediction and annotate as reference spectrum.

    Returns ``(LoadedSpectrum, annotated_df)`` tuple or ``None``.
    """
    try:
        from src.koina_prediction import predict_koina_spectrum, KOINA_AVAILABLE
        if not KOINA_AVAILABLE:
            return None

        import numpy as np
        from src.spectrum_loader import LoadedSpectrum
        from src.peak_annotation import prepare_koina_spectrum_for_plot

        pred_df = predict_koina_spectrum(
            peptidoform_str=peptidoform,
            charge=charge,
            model_id=koina_model,
            collision_energy=koina_ce,
            instrument_type=koina_instrument,
        )
        if pred_df is None or pred_df.empty:
            return None

        mz_arr = pred_df["mz"].values.astype(np.float64)
        int_arr = pred_df["intensity"].values.astype(np.float64)

        # Extract base sequence
        try:
            from psm_utils import Peptidoform as _PF
            base_seq = _PF(peptidoform.split("/")[0]).sequence
        except Exception:
            base_seq = peptidoform.split("/")[0]

        ref_spec = LoadedSpectrum(
            spectrum_id="koina_prediction",
            mz=mz_arr,
            intensity=int_arr,
            precursor_mz=0.0,
            precursor_charge=charge,
            peptidoform=peptidoform,
            base_sequence=base_seq,
            rt=0.0,
            source_file=f"Koina:{koina_model}",
            is_reference=True,
        )

        # Use Koina's native AI-predicted annotations directly
        ref_df = prepare_koina_spectrum_for_plot(pred_df, **annotation_kwargs)
        return (ref_spec, ref_df)
    except Exception as e:
        logger.warning(f"Koina prediction failed for '{peptidoform}': {e}")
        return None


# ============================================================================
# Comparison row builder
# ============================================================================

def build_comparison_row(psm_id, result, query_df, ref_spec, ref_df, metrics):
    """Build a single comparison row dict from query/reference pair.

    *result* is a ``QueryResult`` from ``load_spectra_streaming``.
    *metrics* is the dict returned by ``compute_all_metrics``.
    """
    from src.spectrum_metrics import METRIC_SHORT_TO_DISPLAY

    row: Dict[str, Any] = {
        "psm_id": psm_id,
        "Peptide": result.base_sequence,
        "Modification": result.modification_key,
        "Charge": result.query_spectrum.precursor_charge,
        "Query": Path(result.query_spectrum.source_file).stem,
        "Reference": Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina",
        "query_spectrum_id": result.query_spectrum.spectrum_id,
        "ref_spectrum_id": ref_spec.spectrum_id,
        "query_rt": result.query_spectrum.rt,
        "ref_rt": ref_spec.rt,
    }
    for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
        row[display_name] = metrics.get(short_key, float("nan"))
    for count_key in ("total_q", "total_r", "ann_q", "ann_r", "shared"):
        row[count_key] = metrics.get(count_key, 0)

    return row


# ============================================================================
# process_query_result — core of the shared on_query_ready callback
# ============================================================================

def process_query_result(result, psm_lookup, annotation_kwargs,
                         reference_store, koina_enabled=False,
                         koina_model="ms2pip_timsTOF2024",
                         koina_ce=35, koina_instrument="TIMSTOF"):
    """Process a single query result: annotate, cache refs, compute metrics.

    Parameters
    ----------
    result : QueryResult
        From ``load_spectra_streaming``.
    psm_lookup : dict
        ``{spectrum_lookup_key: psm_id}``
    annotation_kwargs : dict
        Passed through to ``prepare_spectrum_for_plot``.
    reference_store : dict
        Mutable ``{variant_key: [(LoadedSpectrum, df)]}`` -- updated in place.
    koina_enabled, koina_model, koina_ce, koina_instrument :
        Koina fallback settings.

    Returns
    -------
    tuple or None
        ``(psm_id, query_spec, query_df, comparison_rows)`` where
        *comparison_rows* is a list of metric dicts.  Returns ``None``
        when the result cannot be processed (missing PSM, annotation failure).
    """
    from src.spectrum_loader import spectrum_lookup_key
    from src.peak_annotation import prepare_spectrum_for_plot
    from src.spectrum_metrics import compute_all_metrics
    from src.parallel_processor import compute_ion_match_statistics

    key = spectrum_lookup_key(result.query_spectrum)
    psm_id = psm_lookup.get(key)
    if psm_id is None:
        return None

    # Annotate query
    try:
        query_df = prepare_spectrum_for_plot(
            result.query_spectrum.mz,
            result.query_spectrum.intensity,
            result.query_spectrum.peptidoform,
            result.query_spectrum.precursor_charge,
            **annotation_kwargs,
        )
    except Exception as e:
        logger.warning(f"Annotation failed psm_id={psm_id}: {e}")
        return None

    # Annotate and cache references for this variant (once)
    vk = result.variant_key
    if vk not in reference_store and result.reference_spectra:
        ref_list = []
        for ref_spec in result.reference_spectra:
            try:
                ref_df = prepare_spectrum_for_plot(
                    ref_spec.mz, ref_spec.intensity,
                    ref_spec.peptidoform, ref_spec.precursor_charge,
                    **annotation_kwargs,
                )
                ref_list.append((ref_spec, ref_df))
            except Exception as e:
                logger.warning(f"Reference annotation failed: {e}")
        reference_store[vk] = ref_list

    # Koina fallback
    refs = reference_store.get(vk, [])
    if not refs and koina_enabled:
        koina_ref = get_koina_prediction(
            result.query_spectrum.peptidoform,
            result.query_spectrum.precursor_charge,
            annotation_kwargs,
            koina_model=koina_model,
            koina_ce=koina_ce,
            koina_instrument=koina_instrument,
        )
        if koina_ref is not None:
            reference_store[vk] = [koina_ref]
            refs = reference_store[vk]

    # Compute comparison metrics against each reference
    comp_rows: List[Dict[str, Any]] = []
    for ref_spec, ref_df in refs:
        try:
            metrics = compute_all_metrics(query_df, ref_df)
            row = build_comparison_row(
                psm_id, result, query_df, ref_spec, ref_df, metrics,
            )

            # Ion match statistics
            try:
                ims_stats = compute_ion_match_statistics(
                    query_df, result.query_spectrum.peptidoform,
                    result.query_spectrum.precursor_charge,
                    result.query_spectrum.precursor_mz,
                )
                row.update(ims_stats)
            except Exception:
                pass

            comp_rows.append(row)
        except Exception as e:
            logger.warning(f"Metric computation failed: {e}")

    return (psm_id, result.query_spectrum, query_df, comp_rows)


# ============================================================================
# Plot task building
# ============================================================================

def build_plot_tasks(annotated_queries, reference_store, variant_groups,
                     normalize, sanitize_fn, fmt_dirs=None, figsize=(12, 8),
                     dpi=300, html_dir="", plot_unpaired=False):
    """Build list of plot task dicts for parallel/sequential generation.

    Parameters
    ----------
    annotated_queries : list
        ``[(psm_id, LoadedSpectrum, annotated_df), ...]``
    reference_store : dict
        ``{variant_key: [(LoadedSpectrum, annotated_df), ...]}``
    variant_groups : dict
        ``{variant_key: VariantGroup}``
    normalize : bool
    sanitize_fn : callable
        ``sanitize_filename(str) -> str``
    fmt_dirs : dict or None
        ``{"png": "/path", "svg": "/path", ...}`` -- only listed formats are exported.
    figsize : tuple
    dpi : int
    html_dir : str
        Path for HTML export; ``""`` to skip.
    plot_unpaired : bool
        If True, generate single-spectrum plots for queries with no reference.

    Returns
    -------
    list[dict]
    """
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text

    if fmt_dirs is None:
        fmt_dirs = {}

    psm_id_to_vk = build_psm_to_variant_mapping(
        variant_groups, annotated_queries,
    )

    tasks: List[dict] = []
    for pid, query_spec, query_df in sorted(annotated_queries, key=lambda x: x[0]):
        vk = psm_id_to_vk.get(pid)
        refs = reference_store.get(vk, []) if vk else []
        folder_name = sanitize_fn(query_spec.peptidoform)

        if refs:
            for ref_spec, ref_df in refs:
                metrics = compute_all_metrics(query_df, ref_df)
                metrics_text = format_metrics_text(metrics)

                query_label, ref_label = format_spectrum_labels(
                    query_spec, ref_spec,
                )
                safe_name = sanitize_fn(
                    f"{get_source_label(query_spec)}_vs_{get_source_label(ref_spec)}"
                )

                task: Dict[str, Any] = {
                    "type": "mirror",
                    "query_df": query_df,
                    "ref_df": ref_df,
                    "title": f"{query_spec.peptidoform} +{query_spec.precursor_charge}",
                    "query_label": query_label,
                    "ref_label": ref_label,
                    "metrics_text": metrics_text,
                    "normalize": normalize,
                    "folder_name": folder_name,
                    "safe_name": safe_name,
                    "figsize": figsize,
                    "dpi": dpi,
                    "html_dir": html_dir,
                }
                task.update({f"{f}_dir": d for f, d in fmt_dirs.items()})
                tasks.append(task)
        elif plot_unpaired:
            safe_name = sanitize_fn(get_source_label(query_spec))
            query_label = f"{get_source_label(query_spec)} (RT: {query_spec.rt:.1f})"
            task = {
                "type": "single",
                "query_df": query_df,
                "ref_df": None,
                "title": f"{query_spec.peptidoform} +{query_spec.precursor_charge}",
                "query_label": query_label,
                "ref_label": "",
                "metrics_text": "",
                "normalize": normalize,
                "folder_name": folder_name,
                "safe_name": safe_name,
                "figsize": figsize,
                "dpi": dpi,
                "html_dir": "",
            }
            task.update({f"{f}_dir": d for f, d in fmt_dirs.items()})
            tasks.append(task)

    return tasks
