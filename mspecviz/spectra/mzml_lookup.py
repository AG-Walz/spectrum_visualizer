"""mzML spectrum-ID indexing with format-aware fast paths.

CRITICAL: mzML ``index=N`` native IDs are not always 0-based — timsTOF/Bruker
emit 1-based. ``open_indexed_mzml_fast`` auto-detects the offset from the first
spectrum and ``find_spectrum_index`` applies it; removing this offset logic
silently corrupts lookups for Bruker files.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyopenms as oms

from mspecviz.spectra.mzml_native_ids import (
    ANY_DIGITS_RE,
    INDEX_TAG,
    NATIVE_ID_NUM_RE,
    SCAN_TAG,
    SPECTRUM_TAG,
    _extract_id_number,
    clean_spectrum_id,
)

logger = logging.getLogger(__name__)


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """Open mzML and build a targeted ID→position map (skip full scan when possible).

    ``index=N``: O(1) via detected offset. ``scan=N`` with ``needed_ids``: targeted
    lookup. Otherwise: full scan. The returned dict carries three sentinel keys
    (``__format__``, ``__n_spectra__``, ``__index_offset__``) used by
    :func:`find_spectrum_index`.
    """
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))

    n_spectra = indexed_exp.getNrSpectra()
    id_to_index: Dict[str, int] = {}

    # Detect format and index offset from first spectrum
    id_format = "unknown"
    index_offset = 0
    if n_spectra > 0:
        first_id = indexed_exp.getSpectrum(0).getNativeID()
        if INDEX_TAG in first_id:
            id_format = "index"
            # Detect 0-based vs 1-based indexing
            first_num = _extract_id_number(first_id)
            if first_num is not None:
                index_offset = first_num  # 0 for 0-based, 1 for 1-based
        elif SCAN_TAG in first_id:
            id_format = "scan"
        elif SPECTRUM_TAG in first_id:
            id_format = "spectrum"

    id_to_index["__format__"] = id_format
    id_to_index["__n_spectra__"] = n_spectra
    id_to_index["__index_offset__"] = index_offset

    # FAST PATH: "index=N" format - apply offset to convert native ID number
    # to array position. E.g., 1-based "index=1" -> array position 0.
    if id_format == "index":
        logger.debug(f"Fast path: index format detected (offset={index_offset}), "
                     f"skipping full scan ({n_spectra} spectra)")
        return indexed_exp, id_to_index

    # FAST PATH: known needed IDs + numeric format -> targeted lookup
    if needed_ids and id_format in ("scan", "spectrum"):
        found = 0
        for sid in needed_ids:
            clean_id = clean_spectrum_id(sid)
            num = _extract_id_number(clean_id)
            if num is None:
                continue

            # Try direct + off-by-one first (scan numbers often near indices)
            candidates = []
            if 0 <= num < n_spectra:
                candidates.append(num)
            if num > 0 and (num - 1) < n_spectra:
                candidates.append(num - 1)

            matched = False
            for candidate_idx in candidates:
                try:
                    check_spec = indexed_exp.getSpectrum(candidate_idx)
                    check_id = check_spec.getNativeID()
                    if clean_id == clean_spectrum_id(check_id):
                        id_to_index[clean_id] = candidate_idx
                        check_num = _extract_id_number(check_id)
                        if check_num is not None:
                            id_to_index[f"_num_{check_num}"] = candidate_idx
                        found += 1
                        matched = True
                        break
                except Exception:
                    continue

            if not matched:
                # Binary search: scan numbers are typically monotonic
                lo, hi = 0, n_spectra - 1
                while lo <= hi:
                    mid = (lo + hi) // 2
                    try:
                        mid_spec = indexed_exp.getSpectrum(mid)
                        mid_id = mid_spec.getNativeID()
                        mid_num = _extract_id_number(mid_id)
                        if mid_num is None:
                            break
                        if mid_num == num:
                            id_to_index[clean_id] = mid
                            id_to_index[f"_num_{mid_num}"] = mid
                            found += 1
                            break
                        elif mid_num < num:
                            lo = mid + 1
                        else:
                            hi = mid - 1
                    except Exception:
                        break

        logger.debug(f"Targeted lookup: found {found}/{len(needed_ids)} spectra in {Path(mzml_path).name}")

        # If we found enough, skip the slow fallback
        if found > 0:
            return indexed_exp, id_to_index

    # SLOW FALLBACK: full scan (unknown format or all targeted lookups missed)
    logger.debug(f"Full scan fallback for {Path(mzml_path).name} ({n_spectra} spectra)")
    for i in range(n_spectra):
        try:
            spec = indexed_exp.getSpectrum(i)
            native_id = spec.getNativeID()
            id_to_index[native_id] = i
            num_match = ANY_DIGITS_RE.search(native_id)
            if num_match:
                id_to_index[f"_num_{num_match.group(1)}"] = i
        except Exception:
            continue

    return indexed_exp, id_to_index


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int,
) -> Optional[int]:
    """Resolve a spectrum ID to its array position via the map from ``open_indexed_mzml_fast``.

    Tries direct lookup, then ``(scan|index|spectrum)=N`` with offset adjustment,
    then a bare-digit fallback.
    """
    clean_id = clean_spectrum_id(spectrum_id)
    if not clean_id:
        return None

    id_format = id_to_index.get("__format__", "unknown")
    index_offset = id_to_index.get("__index_offset__", 0)

    # Strategy 1: direct lookup
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    # Strategy 2: structured ID with (scan|index|spectrum)=N
    scan_match = NATIVE_ID_NUM_RE.search(clean_id)
    if scan_match:
        num = int(scan_match.group(1))

        if id_format == "index":
            arr_pos = num - index_offset
            if 0 <= arr_pos < n_spectra:
                return arr_pos

        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

    # Strategy 3: bare number fallback
    num_match = ANY_DIGITS_RE.search(clean_id)
    if num_match:
        num = int(num_match.group(1))
        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

        if id_format == "index":
            arr_pos = num - index_offset
            if 0 <= arr_pos < n_spectra:
                return arr_pos

    return None
