"""Per-query annotation phases used by ``process_query_result``.

Each phase is a small helper: annotate the query, cache references, attempt a
Koina fallback, build metric rows. Extracted so ``query_processing`` stays
readable and each phase can be unit-tested independently.
"""
from typing import Any

import numpy as np
from psm_utils import Peptidoform

from mspecviz.annotation.comparison_metrics import METRIC_SHORT_TO_DISPLAY
from mspecviz.annotation.peak_matching import (
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from mspecviz.common.defaults import DEFAULTS
from mspecviz.spectra.koina_prediction import KOINA_AVAILABLE, predict_koina_spectrum
from mspecviz.spectra.loaded_spectrum import LoadedSpectrum, source_label


def get_koina_prediction(
    peptidoform, charge, annotation_kwargs,
    koina_model: str = DEFAULTS.koina_model,
    koina_collision_energy: int = DEFAULTS.koina_collision_energy,
    koina_instrument_type: str = DEFAULTS.koina_instrument_type,
):
    """Get Koina AI prediction and annotate as reference spectrum. """
    # Availability guard: bail early if Koina extra isn't installed
    if not KOINA_AVAILABLE:
        return None

    # Predict: call Koina with model + CE + instrument, bail on empty frame
    predicted_frame = predict_koina_spectrum(
        peptidoform_str=peptidoform,
        charge=charge,
        model_id=koina_model,
        collision_energy=koina_collision_energy,
        instrument_type=koina_instrument_type,
    )
    if predicted_frame is None or predicted_frame.empty:
        return None

    # Wrap prediction: build a synthetic LoadedSpectrum tagged as reference
    reference_spectrum = LoadedSpectrum(
        spectrum_id="koina_prediction",
        mz=predicted_frame["mz"].to_numpy(dtype=np.float64),
        intensity=predicted_frame["intensity"].to_numpy(dtype=np.float64),
        precursor_mz=0.0,
        precursor_charge=charge,
        peptidoform=peptidoform,
        base_sequence=Peptidoform(peptidoform).sequence,
        rt=0.0,
        source_file=f"Koina:{koina_model}",
        is_reference=True,
    )

    # Annotate: reuse Koina's native AI-predicted ion labels directly
    reference_frame = prepare_koina_spectrum_for_plot(predicted_frame, **annotation_kwargs)
    return (reference_spectrum, reference_frame)


def _annotate_query_spectrum(result, annotation_kwargs):
    """Annotate the query spectrum (peaks + theoretical fragments)."""
    return prepare_spectrum_for_plot(
        result.query_spectrum.mz,
        result.query_spectrum.intensity,
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        **annotation_kwargs,
    )


def _cache_reference_spectra(result, reference_store, annotation_kwargs):
    """Annotate and cache reference spectra for this variant (once)."""
    variant_key = result.variant_key
    if variant_key in reference_store or not result.reference_spectra:
        return
    ref_list = []
    for reference_spectrum in result.reference_spectra:
        reference_frame = prepare_spectrum_for_plot(
            reference_spectrum.mz, reference_spectrum.intensity,
            reference_spectrum.peptidoform, reference_spectrum.precursor_charge,
            **annotation_kwargs,
        )
        ref_list.append((reference_spectrum, reference_frame))
    reference_store[variant_key] = ref_list


def _try_koina_fallback(result, reference_store, annotation_kwargs,
                        koina_model, koina_collision_energy, koina_instrument_type):
    """If no experimental references, try a Koina AI prediction. Returns the ref list."""
    variant_key = result.variant_key
    references = reference_store.get(variant_key, [])
    if references:
        return references
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
    return references


def build_comparison_row(psm_id, result, reference_spectrum, reference_frame, metrics):
    """Build a single comparison row dict from query/reference pair.

    *result* is a ``QueryResult`` from ``load_spectra_streaming``.
    *metrics* is the dict returned by ``compute_all_metrics``.
    """
    # Identity columns: per-PSM metadata + query/reference labels and RTs
    row: dict[str, Any] = {
        "psm_id": psm_id,
        "Peptide": result.query_spectrum.base_sequence,
        "Modification": result.variant_key if result.variant_key != result.query_spectrum.base_sequence else "",
        "Charge": result.query_spectrum.precursor_charge,
        "Query": source_label(result.query_spectrum),
        "Reference": source_label(reference_spectrum),
        "query_spectrum_id": result.query_spectrum.spectrum_id,
        "ref_spectrum_id": reference_spectrum.spectrum_id,
        "query_rt": result.query_spectrum.rt,
        "ref_rt": reference_spectrum.rt,
    }
    # Metric columns: short keys -> display names, NaN for missing values
    for short_key, display_name in METRIC_SHORT_TO_DISPLAY.items():
        row[display_name] = metrics.get(short_key, float("nan"))
    # Count columns: peak/annotation tallies on both sides plus shared
    for count_key in (
        "query_peak_count", "reference_peak_count",
        "query_annotation_count", "reference_annotation_count", "shared",
    ):
        row[count_key] = metrics.get(count_key, 0)

    return row
