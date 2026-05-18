"""Build the spectrum-ID -> array-position index for one mzML (used by ``mzml_workers``).

CRITICAL: ``index=N`` is not always 0-based (Bruker timsTOF emits 1-based) — offset auto-detected."""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyopenms as oms

from core.spectra.mzml_native_ids import (
    INDEX_TAG,
    SCAN_TAG,
    SPECTRUM_TAG,
    clean_spectrum_id,
    extract_id_number,
)

logger = logging.getLogger(__name__)


def _detect_id_format(indexed_exp: oms.OnDiscMSExperiment, n_spectra: int) -> Tuple[str, int]:
    """Probe the first spectrum's native ID; return ``(format, offset)`` where ``offset`` is 0 or 1 for ``index=N`` only (Bruker timsTOF gotcha)."""
    if n_spectra == 0:
        return "unknown", 0
    first_id = indexed_exp.getSpectrum(0).getNativeID()
    if INDEX_TAG in first_id:
        first_number = extract_id_number(first_id)
        return "index", (first_number if first_number is not None else 0)
    if SCAN_TAG in first_id:
        return "scan", 0
    if SPECTRUM_TAG in first_id:
        return "spectrum", 0
    return "unknown", 0


def _try_candidate_indices(
    indexed_exp: oms.OnDiscMSExperiment, clean_id: str, candidates: List[int],
    id_to_index: Dict[str, int],
) -> bool:
    """Check a few likely array positions for a matching native ID; record the first hit and return ``True``."""
    for candidate_index in candidates:
        check_id = indexed_exp.getSpectrum(candidate_index).getNativeID()
        if clean_id != clean_spectrum_id(check_id):
            continue
        id_to_index[clean_id] = candidate_index
        check_number = extract_id_number(check_id)
        if check_number is not None:
            id_to_index[f"_num_{check_number}"] = candidate_index
        return True
    return False


def _binary_search_for_number(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, target_number: int,
    clean_id: str, id_to_index: Dict[str, int],
) -> bool:
    """Binary-search for a scan number across spectra (scan numbers are typically monotonic)."""
    low, high = 0, n_spectra - 1
    while low <= high:
        mid = (low + high) // 2
        try:
            mid_id = indexed_exp.getSpectrum(mid).getNativeID()
        except Exception:
            return False
        mid_number = extract_id_number(mid_id)
        if mid_number is None:
            return False
        if mid_number == target_number:
            id_to_index[clean_id] = mid
            id_to_index[f"_num_{mid_number}"] = mid
            return True
        if mid_number < target_number:
            low = mid + 1
        else:
            high = mid - 1
    return False


def _resolve_targeted_id(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, spectrum_id: str,
    id_to_index: Dict[str, int],
) -> bool:
    """Locate one needed spectrum ID — try off-by-one candidates first, then binary search by scan number."""
    clean_id = clean_spectrum_id(spectrum_id)
    extracted_number = extract_id_number(clean_id)
    if extracted_number is None:
        return False
    candidates = [c for c in (extracted_number, extracted_number - 1) if 0 <= c < n_spectra]
    if _try_candidate_indices(indexed_exp, clean_id, candidates, id_to_index):
        return True
    return _binary_search_for_number(
        indexed_exp, n_spectra, extracted_number, clean_id, id_to_index,
    )


def _build_full_scan_index(
    indexed_exp: oms.OnDiscMSExperiment, n_spectra: int, id_to_index: Dict[str, int],
) -> None:
    """Last-resort: walk every spectrum and index it by native ID."""
    for index in range(n_spectra):
        try:
            native_id = indexed_exp.getSpectrum(index).getNativeID()
        except Exception:
            continue
        id_to_index[native_id] = index
        number = extract_id_number(native_id)
        if number is not None:
            id_to_index[f"_num_{number}"] = index


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """Open one mzML and build its ID -> array-position map, skipping a full scan when possible.

    Returns ``(experiment, id_to_index)``; dict carries ``__format__``/``__n_spectra__``/``__index_offset__`` sentinels."""
    # Open & detect: probe the first spectrum's native-ID format and offset
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))
    n_spectra = indexed_exp.getNrSpectra()

    id_format, index_offset = _detect_id_format(indexed_exp, n_spectra)
    id_to_index: Dict[str, int] = {
        "__format__": id_format,
        "__n_spectra__": n_spectra,
        "__index_offset__": index_offset,
    }

    # Fast path: index=N maps to array position via the detected offset
    if id_format == "index":
        logger.debug(f"Fast path: index format detected (offset={index_offset}), "
                     f"skipping full scan ({n_spectra} spectra)")
        return indexed_exp, id_to_index

    # Targeted lookup: resolve only the needed scan/spectrum IDs when supplied
    if needed_ids and id_format in ("scan", "spectrum"):
        found = sum(
            _resolve_targeted_id(indexed_exp, n_spectra, sid, id_to_index)
            for sid in needed_ids
        )
        logger.debug(
            f"Targeted lookup: found {found}/{len(needed_ids)} spectra in "
            f"{Path(mzml_path).name}"
        )
        if found > 0:
            return indexed_exp, id_to_index

    # Fallback: walk every spectrum and index by native ID
    logger.debug(f"Full scan fallback for {Path(mzml_path).name} ({n_spectra} spectra)")
    _build_full_scan_index(indexed_exp, n_spectra, id_to_index)
    return indexed_exp, id_to_index


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int,
) -> Optional[int]:
    """Resolve one spectrum ID to its array position using the map from ``open_indexed_mzml_fast``; ``None`` if unknown."""
    clean_id = clean_spectrum_id(spectrum_id)
    if not clean_id:
        return None
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    extracted = extract_id_number(clean_id)
    if extracted is None:
        return None

    sentinel = id_to_index.get(f"_num_{extracted}")
    if sentinel is not None:
        return sentinel

    if id_to_index.get("__format__") == "index":
        offset = id_to_index.get("__index_offset__", 0)
        position = extracted - offset
        if 0 <= position < n_spectra:
            return position

    return None
