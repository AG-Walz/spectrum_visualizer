"""
mzML spectrum-ID indexing with format-aware fast paths.

Extracted from ``spectrum_loader`` so the lookup logic (which has non-obvious
edge cases around 1-based Bruker/timsTOF ``index=N`` native IDs) lives in a
focused module that can be unit-tested independently.

Key concern: mzML ``index=N`` native IDs are **not necessarily 0-based**.
timsTOF/Bruker converters emit 1-based IDs, so ``index=1`` maps to array
position 0. ``open_indexed_mzml_fast`` auto-detects the offset from the first
spectrum; ``find_spectrum_index`` applies it on lookup. Do not remove or
reorder this logic without reading ``docs/`` and running integration tests.
"""
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyopenms as oms

logger = logging.getLogger(__name__)


def _extract_id_number(native_id: str) -> Optional[int]:
    """Extract the numeric value from a native ID string."""
    m = re.search(r'(?:scan|index|spectrum)=(\d+)', native_id)
    if m:
        return int(m.group(1))
    m = re.search(r'(\d+)', native_id)
    if m:
        return int(m.group(1))
    return None


def clean_spectrum_id(spectrum_id: str) -> str:
    """Clean a spectrum ID for lookup (strips whitespace, extracts canonical form)."""
    if not spectrum_id:
        return ""

    if isinstance(spectrum_id, bytes):
        spectrum_id = spectrum_id.decode()

    clean_id = spectrum_id.strip()

    # If it contains "controllerType" (Thermo format), keep it
    if "controllerType" in clean_id:
        return clean_id

    # For simple scan=N format, ensure clean extraction
    for tag in ["scan=", "index=", "spectrum="]:
        if tag in clean_id:
            idx = clean_id.find(tag)
            rest = clean_id[idx:]
            parts = rest.split()
            return parts[0]

    return clean_id


def open_indexed_mzml_fast(
    mzml_path: str,
    needed_ids: Optional[List[str]] = None
) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """
    Open mzML with targeted (not full) index building.

    Optimizations over a naive full-scan:
    - For ``index=N`` format: O(1) array lookup after offset detection.
    - For ``scan=N`` with known *needed_ids*: targeted lookup + binary search.
    - Falls back to full scan when format or *needed_ids* are unavailable.

    Args:
        mzml_path: Path to mzML file.
        needed_ids: Optional list of spectrum IDs we need. If provided, enables
            targeted lookup instead of a full scan.

    Returns:
        ``(OnDiscMSExperiment, id_to_index)`` where ``id_to_index`` also carries
        three special keys: ``__format__``, ``__n_spectra__``, ``__index_offset__``.
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
        if "index=" in first_id:
            id_format = "index"
            # Detect 0-based vs 1-based indexing
            first_num = _extract_id_number(first_id)
            if first_num is not None:
                index_offset = first_num  # 0 for 0-based, 1 for 1-based
        elif "scan=" in first_id:
            id_format = "scan"
        elif "spectrum=" in first_id:
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
            num_match = re.search(r'(\d+)', native_id)
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
    """Find spectrum array index using the lookup produced by ``open_indexed_mzml_fast``.

    Strategies, in order:
        1. Direct dict lookup on the cleaned ID.
        2. Extract numeric value; for ``index=`` format apply offset; otherwise
           check the ``_num_N`` cache; fall back to offset-adjusted array access.
        3. Last-ditch: extract any digits from the ID and try the same paths.
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
    scan_match = re.search(r'(?:scan|index|spectrum)=(\d+)', clean_id)
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
    num_match = re.search(r'(\d+)', clean_id)
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
