"""Streaming mzML loader — picklable worker.

Only holds the module-level worker that can be shipped across the
ProcessPoolExecutor boundary. The orchestration logic (map-building,
phase ordering, callbacks) lives in ``core.spectra.spectrum_streaming``.
"""
import logging
from typing import List, Tuple

from psm_utils import PSM

from core.spectra.loaded_spectrum import LoadResult, LoadedSpectrum
from core.spectra.mzml_lookup import find_spectrum_index, open_indexed_mzml_fast

logger = logging.getLogger(__name__)


def _build_loaded_spectrum(
    spectrum, psm: PSM, display_name: str, is_reference: bool,
) -> LoadedSpectrum:
    """Wrap a pyopenms ``MSSpectrum`` in our ``LoadedSpectrum`` dataclass."""
    mz_array, intensity_array = spectrum.get_peaks()
    psm_charge = psm.get_precursor_charge() or 0
    # Precursor: pull m/z and charge from the spectrum, falling back to psm_charge
    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        precursor_mz = prec.getMZ()
        precursor_charge = prec.getCharge() if prec.getCharge() > 0 else psm_charge
    else:
        precursor_mz = 0.0
        precursor_charge = psm_charge
    return LoadedSpectrum(
        spectrum_id=psm.spectrum_id,
        mz=mz_array,
        intensity=intensity_array,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        peptidoform=psm.peptidoform.modified_sequence,
        base_sequence=psm.peptidoform.sequence.upper(),
        source_file=display_name,
        is_reference=is_reference,
        rt=spectrum.getRT(),
        match=psm,
    )


def _load_spectra_from_mzml(
    mzml_path_str: str,
    display_name: str,
    matches_for_file: List[Tuple[PSM, str]],
    is_reference: bool,
) -> LoadResult:
    """Worker: open one mzML and extract matched spectra. Module-level for pickling."""
    # Setup: per-file accumulators for loaded spectra and errors
    loaded: List[Tuple[str, LoadedSpectrum]] = []
    errors: List[str] = []

    try:
        # Open & index: build targeted ID→position map for the needed spectra
        needed_ids = [psm.spectrum_id for psm, _ in matches_for_file]
        indexed_exp, id_to_index = open_indexed_mzml_fast(mzml_path_str, needed_ids=needed_ids)
        n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

        # Resolve & wrap: look up each PSM's spectrum and build LoadedSpectrum
        for psm, variant_key in matches_for_file:
            try:
                spectrum_index = find_spectrum_index(psm.spectrum_id, id_to_index, n_spectra)
                if spectrum_index is None:
                    errors.append(f"Spectrum not found: {psm.spectrum_id} in {display_name}")
                    continue
                spectrum = indexed_exp.getSpectrum(spectrum_index)
                loaded.append((variant_key, _build_loaded_spectrum(
                    spectrum, psm, display_name, is_reference,
                )))
            except Exception as error:
                errors.append(f"Error loading {psm.spectrum_id}: {error}")

    except Exception as error:
        errors.append(f"Error opening {display_name}: {error}")

    # Package: counts + errors into LoadResult for the orchestrator
    return LoadResult(
        display_name=display_name,
        is_reference=is_reference,
        loaded=loaded,
        loaded_count=len(loaded),
        not_found_count=len(matches_for_file) - len(loaded),
        errors=errors,
    )
