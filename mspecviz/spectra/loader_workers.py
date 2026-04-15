"""Streaming mzML loader — picklable worker.

Only holds the module-level worker that can be shipped across the
ProcessPoolExecutor boundary. The orchestration logic (map-building,
phase ordering, callbacks) lives in ``mspecviz.spectra.loader_orchestration``.
"""
import logging
from typing import List, Tuple

from mspecviz.search.models import PeptideMatch
from mspecviz.spectra.loader_models import LoadResult, LoadedSpectrum
from mspecviz.spectra.mzml_lookup import find_spectrum_index, open_indexed_mzml_fast

logger = logging.getLogger(__name__)


def _load_spectra_from_mzml(
    mzml_path_str: str,
    display_name: str,
    matches_for_file: List[Tuple["PeptideMatch", str]],
    is_reference: bool,
) -> LoadResult:
    """Worker: open one mzML and extract matched spectra. Module-level for pickling."""
    loaded: List[Tuple[str, LoadedSpectrum]] = []
    errors: List[str] = []
    loaded_count = 0
    not_found_count = 0

    try:
        needed_ids = [match.spectrum_id for match, _ in matches_for_file]
        indexed_exp, id_to_index = open_indexed_mzml_fast(mzml_path_str, needed_ids=needed_ids)
        n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

        for match, variant_key in matches_for_file:
            try:
                idx = find_spectrum_index(match.spectrum_id, id_to_index, n_spectra)
                if idx is None:
                    not_found_count += 1
                    errors.append(f"Spectrum not found: {match.spectrum_id} in {display_name}")
                    continue

                spectrum = indexed_exp.getSpectrum(idx)
                mz_array, intensity_array = spectrum.get_peaks()

                precursor_mz = 0.0
                precursor_charge = match.charge
                if spectrum.getPrecursors():
                    prec = spectrum.getPrecursors()[0]
                    precursor_mz = prec.getMZ()
                    if prec.getCharge() > 0:
                        precursor_charge = prec.getCharge()

                loaded.append((variant_key, LoadedSpectrum(
                    spectrum_id=match.spectrum_id,
                    mz=mz_array,
                    intensity=intensity_array,
                    precursor_mz=precursor_mz,
                    precursor_charge=precursor_charge,
                    peptidoform=match.peptidoform,
                    base_sequence=match.base_sequence,
                    source_file=display_name,
                    is_reference=is_reference,
                    rt=spectrum.getRT(),
                    match=match,
                )))
                loaded_count += 1

            except Exception as e:
                not_found_count += 1
                errors.append(f"Error loading {match.spectrum_id}: {e}")

    except Exception as e:
        errors.append(f"Error opening {display_name}: {e}")

    return LoadResult(
        display_name=display_name,
        is_reference=is_reference,
        loaded=loaded,
        loaded_count=loaded_count,
        not_found_count=not_found_count,
        errors=errors,
    )
