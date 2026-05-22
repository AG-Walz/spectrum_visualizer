"""Per-file mzML loading worker — opens one mzML and extracts the spectra named by a batch of PSMs.
This module holds the unit of work farmed out by the spectrum-streaming
stage: each call to _load_spectra_from_mzml handles exactly one mzML file
and the subset of PSMs whose spectra live in it. The result is a
LoadResult bundling the wrapped LoadedSpectrum objects, per-spectrum
counts, and any errors encountered.
Everything here must stay picklable and free of process-global state — the
worker runs in a separate process under ProcessPoolExecutor and only plain
Python / numpy values may cross the process boundary back to the caller.
(called by spectra.spectrum_streaming via helpers.run_tasks; uses
spectra.mzml_lookup to open and index the file, and
spectra.loaded_spectrum for the LoadedSpectrum / LoadResult result types.)
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
    """Wrap a pyopenms MSSpectrum and its matched PSM into a single LoadedSpectrum.
    Pulls the peak arrays, precursor m/z and charge, retention time, and
    peptidoform off the raw objects. When the spectrum carries no precursor
    block — or its charge is unset — the precursor charge falls back to the
    PSM's charge, and precursor_mz defaults to 0.0. Returns the fully
    populated LoadedSpectrum; runs inside a worker process, so the result
    must stay picklable.
    """
    mz_array, intensity_array = spectrum.get_peaks()
    # PSM charge is the fallback precursor charge when the spectrum lacks one
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
    """Worker entry point — open one mzML and load every spectrum named in matches_for_file.
    matches_for_file is a list of (PSM, variant_key) pairs; each PSM's
    spectrum_id is resolved to an array position via a targeted index, the
    spectrum is fetched, and it is wrapped into a LoadedSpectrum tagged with
    its variant key. Errors are collected per spectrum (and once for the whole
    file if it fails to open) rather than raised — a failure on one spectrum or
    file never aborts the batch.
    Returns a LoadResult carrying loaded as a list of
    (variant_key, LoadedSpectrum) pairs, plus the loaded / not-found counts
    and the error list. Runs inside a ProcessPoolExecutor worker, so the
    return value must remain picklable.
    """
    # Per-file accumulators — successfully loaded spectra and error strings
    loaded: List[Tuple[str, LoadedSpectrum]] = []
    errors: List[str] = []

    try:
        # Open the file once and build an ID -> array-position map limited to
        # the spectra this batch actually needs — avoids a full-file scan.
        needed_ids = [psm.spectrum_id for psm, _ in matches_for_file]
        indexed_exp, id_to_index = open_indexed_mzml_fast(mzml_path_str, needed_ids=needed_ids)
        n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())

        # Resolve each PSM's spectrum ID to a position, fetch it, and wrap it
        for psm, variant_key in matches_for_file:
            try:
                spectrum_index = find_spectrum_index(psm.spectrum_id, id_to_index, n_spectra)
                # Missing ID: record and skip — one absent spectrum is not fatal
                if spectrum_index is None:
                    errors.append(f"Spectrum not found: {psm.spectrum_id} in {display_name}")
                    continue
                spectrum = indexed_exp.getSpectrum(spectrum_index)
                loaded.append((variant_key, _build_loaded_spectrum(
                    spectrum, psm, display_name, is_reference,
                )))
            # Per-spectrum failure: keep going so the rest of the batch loads
            except Exception as error:
                errors.append(f"Error loading {psm.spectrum_id}: {error}")

    # Whole-file failure (e.g. unreadable mzML): record and return an empty result
    except Exception as error:
        errors.append(f"Error opening {display_name}: {error}")

    # Bundle counts and errors into a LoadResult for the streaming orchestrator
    return LoadResult(
        display_name=display_name,
        is_reference=is_reference,
        loaded=loaded,
        loaded_count=len(loaded),
        not_found_count=len(matches_for_file) - len(loaded),
        errors=errors,
    )
