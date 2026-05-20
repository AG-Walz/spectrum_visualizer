"""Per-query annotation phases — annotate query, cache references, fall back to Koina (called by ``pipeline.query_processing.process_query_result``)."""
import numpy as np
from psm_utils import Peptidoform

from core.annotation.peak_matching import (
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.defaults import DEFAULTS
from core.spectra.koina_prediction import predict_koina_spectrum
from core.spectra.loaded_spectrum import LoadedSpectrum


def get_koina_prediction(
    peptidoform, charge, annotation_kwargs,
    koina_model: str = DEFAULTS.koina_model,
    koina_collision_energy: int = DEFAULTS.koina_collision_energy,
    koina_instrument_type: str = DEFAULTS.koina_instrument_type,
):
    """Predict a spectrum with Koina and wrap it as an annotated reference; return ``(spectrum, frame)`` or ``None``."""
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
    """Match the query spectrum's peaks against theoretical fragment ions; return the annotated frame."""
    return prepare_spectrum_for_plot(
        result.query_spectrum.mz,
        result.query_spectrum.intensity,
        result.query_spectrum.peptidoform,
        result.query_spectrum.precursor_charge,
        precursor_mz=result.query_spectrum.precursor_mz,
        **annotation_kwargs,
    )


def _cache_reference_spectra(result, reference_store, annotation_kwargs):
    """Annotate every reference spectrum for this peptide variant and stash them in ``reference_store`` (once per variant)."""
    variant_key = result.variant_key
    if variant_key in reference_store or not result.reference_spectra:
        return
    ref_list = []
    for reference_spectrum in result.reference_spectra:
        reference_frame = prepare_spectrum_for_plot(
            reference_spectrum.mz, reference_spectrum.intensity,
            reference_spectrum.peptidoform, reference_spectrum.precursor_charge,
            precursor_mz=reference_spectrum.precursor_mz,
            **annotation_kwargs,
        )
        ref_list.append((reference_spectrum, reference_frame))
    reference_store[variant_key] = ref_list


def _try_koina_fallback(result, reference_store, annotation_kwargs,
                        koina_model, koina_collision_energy, koina_instrument_type):
    """When no experimental reference exists for this variant, request a Koina prediction; return the resolved list."""
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
