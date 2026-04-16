"""Per-query annotation phases used by ``process_query_result``.

Each phase is a small helper: annotate the query, cache references, attempt a
Koina fallback, build metric rows. Extracted so ``orchestration.py`` stays
readable and each phase can be unit-tested independently.
"""
import logging
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from psm_utils import Peptidoform

from mspecviz.annotation.metrics_set import METRIC_SHORT_TO_DISPLAY, compute_all_metrics
from mspecviz.annotation.peak_matching import (
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from mspecviz.spectra.koina_prediction import KOINA_AVAILABLE, predict_koina_spectrum
from mspecviz.spectra.loader_models import LoadedSpectrum

logger = logging.getLogger("mspecviz.pipeline")


def get_koina_prediction(peptidoform, charge, annotation_kwargs,
                         koina_model="ms2pip_timsTOF2024",
                         koina_collision_energy=35, koina_instrument_type="TIMSTOF"):
    """Get Koina AI prediction and annotate as reference spectrum.

    Returns ``(LoadedSpectrum, annotated_df)`` tuple or ``None``.
    """
    if not KOINA_AVAILABLE:
        return None

    pred_df = predict_koina_spectrum(
        peptidoform_str=peptidoform,
        charge=charge,
        model_id=koina_model,
        collision_energy=koina_collision_energy,
        instrument_type=koina_instrument_type,
    )
    if pred_df is None or pred_df.empty:
        return None

    mz_arr = pred_df["mz"].values.astype(np.float64)
    int_arr = pred_df["intensity"].values.astype(np.float64)

    base_seq = Peptidoform(peptidoform.split("/")[0]).sequence

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


def _annotate_query_spectrum(result, annotation_kwargs, psm_id):
    """Phase 1 of ``process_query_result``: annotate the query spectrum."""
    try:
        return prepare_spectrum_for_plot(
            result.query_spectrum.mz,
            result.query_spectrum.intensity,
            result.query_spectrum.peptidoform,
            result.query_spectrum.precursor_charge,
            **annotation_kwargs,
        )
    except Exception as e:
        logger.warning(f"Annotation failed psm_id={psm_id}: {e}")
        return None


def _cache_reference_spectra(result, reference_store, annotation_kwargs):
    """Phase 2: annotate and cache references for this variant (once)."""
    variant_key = result.variant_key
    if variant_key in reference_store or not result.reference_spectra:
        return
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
    reference_store[variant_key] = ref_list


def _try_koina_fallback(result, reference_store, annotation_kwargs,
                        koina_model, koina_collision_energy, koina_instrument_type):
    """Phase 3: if no experimental refs, try a Koina AI prediction. Returns the ref list."""
    variant_key = result.variant_key
    refs = reference_store.get(variant_key, [])
    if refs:
        return refs
    koina_ref = get_koina_prediction(
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        annotation_kwargs,
        koina_model=koina_model,
        koina_collision_energy=koina_collision_energy,
        koina_instrument_type=koina_instrument_type,
    )
    if koina_ref is not None:
        reference_store[variant_key] = [koina_ref]
        return reference_store[variant_key]
    return refs


def build_comparison_row(psm_id, result, query_df, ref_spec, ref_df, metrics):
    """Build a single comparison row dict from query/reference pair.

    *result* is a ``QueryResult`` from ``load_spectra_streaming``.
    *metrics* is the dict returned by ``compute_all_metrics``.
    """
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


def _build_metric_rows(psm_id, result, query_df, refs):
    """Phase 4: compute metrics + build summary rows for each query/ref pair."""
    rows: List[Dict[str, Any]] = []
    for ref_spec, ref_df in refs:
        try:
            metrics = compute_all_metrics(query_df, ref_df)
            rows.append(build_comparison_row(
                psm_id, result, query_df, ref_spec, ref_df, metrics,
            ))
        except Exception as e:
            logger.warning(f"Metric computation failed: {e}")
    return rows
