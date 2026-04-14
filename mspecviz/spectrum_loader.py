"""
Spectrum loading from mzML files.

Streams spectra with reference-first ordering so downstream consumers can
compute query/reference metrics as each query arrives. Uses the pre-resolved
mzML mappings from the search phase. All indexing logic lives in
``mspecviz.mzml_indexing``; parallel/sequential dispatch lives in
``mspecviz.common.parallel``.
"""
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from mspecviz.common.parallel import run_tasks
from mspecviz.mzml_indexing import find_spectrum_index, open_indexed_mzml_fast
from mspecviz.path_utils import disambiguate_paths
from mspecviz.peptide_search import PeptideMatch, VariantGroup

logger = logging.getLogger(__name__)


@dataclass
class LoadedSpectrum:
    """A loaded spectrum with metadata."""
    spectrum_id: str
    mz: np.ndarray
    intensity: np.ndarray
    precursor_mz: float
    precursor_charge: int
    peptidoform: str
    base_sequence: str
    source_file: str
    is_reference: bool = False
    rt: float = 0.0
    match: Optional[PeptideMatch] = None


@dataclass
class StreamingResult:
    """Result from streaming spectrum loader - one query spectrum with its matching references."""
    query_spectrum: LoadedSpectrum
    reference_spectra: List[LoadedSpectrum]
    variant_key: str
    base_sequence: str
    modification_key: str
    query_index: int
    total_queries: int


@dataclass
class LoadResult:
    """Per-file output from ``_load_spectra_from_mzml``.

    Replaces an anonymous 6-tuple. *loaded* is a list of
    ``(variant_key, LoadedSpectrum)`` pairs so the caller can bucket spectra
    by variant without re-scanning.
    """
    display_name: str
    is_reference: bool
    loaded: List[Tuple[str, LoadedSpectrum]]
    loaded_count: int
    not_found_count: int
    errors: List[str]


def build_mzml_to_matches_map(
    variant_groups: Dict[str, VariantGroup],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PeptideMatch, str]]]:
    """Map each mzML path key to the list of (match, variant_key) pairs to load.

    Uses ``mzml_resolved_path`` when available (handles same-name files from
    different directories). Falls back to bare-filename matching; ambiguous
    matches pick the first candidate and log a warning.
    """
    # Bare-name fallback: lowercase name -> list of path keys
    name_to_keys: Dict[str, List[str]] = {}
    for key, path in available_mzml_paths.items():
        name_to_keys.setdefault(path.name.lower(), []).append(key)

    mzml_to_matches: Dict[str, List[Tuple[PeptideMatch, str]]] = {}

    for variant_key, variant_group in variant_groups.items():
        for match in variant_group.matches:
            if not match.mzml_filename:
                continue

            # Prefer resolved path (handles same-name files)
            path_key = None
            if match.mzml_resolved_path:
                resolved_key = str(Path(match.mzml_resolved_path).resolve()).lower()
                if resolved_key in available_mzml_paths:
                    path_key = resolved_key

            # Fallback: bare filename match
            if path_key is None:
                candidates = name_to_keys.get(match.mzml_filename.lower(), [])
                if candidates:
                    path_key = candidates[0]
                    if len(candidates) > 1:
                        paths = [str(available_mzml_paths[c]) for c in candidates]
                        logger.warning(
                            f"Ambiguous mzML match for '{match.mzml_filename}': "
                            f"{len(candidates)} files with same name. "
                            f"Using '{paths[0]}'. All candidates: {paths}"
                        )

            if path_key:
                mzml_to_matches.setdefault(path_key, []).append((match, variant_key))

    return mzml_to_matches


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


def _build_file_tasks(file_keys, mzml_to_matches, mzml_path_map, mzml_disambig):
    """Build (path_str, display_name, matches) tuples for files that have matches."""
    tasks = []
    for path_key in file_keys:
        if path_key not in mzml_to_matches:
            continue
        mzml_path = mzml_path_map[path_key]
        display_name = mzml_disambig.get(str(mzml_path), mzml_path.name)
        tasks.append((str(mzml_path), display_name, mzml_to_matches[path_key]))
    return tasks


def _execute_load_tasks(file_tasks, is_reference, workers,
                        on_file_start, on_file_complete, error_callback) -> List[LoadResult]:
    """Load spectra from mzML files, sequentially or in parallel."""
    kind = "reference" if is_reference else "query"
    if not file_tasks:
        return []

    if workers > 1 and len(file_tasks) > 1:
        logger.info(f"Parallel {kind} loading: {len(file_tasks)} files, {workers} workers")

    # Tack the fixed is_reference onto each task tuple
    tasks = [(p, d, m, is_reference) for p, d, m in file_tasks]
    all_results: List[LoadResult] = []

    def _on_start(task):
        if on_file_start:
            on_file_start(task[1], is_reference)  # task[1] = display_name

    def _on_result(task, result: LoadResult):
        if error_callback:
            for err in result.errors:
                error_callback(err)
        all_results.append(result)
        if on_file_complete:
            on_file_complete(result.display_name, result.loaded_count, result.not_found_count)

    def _on_error(task, exc):
        display_name = task[1]
        if error_callback:
            error_callback(f"Worker crashed loading {kind} {display_name}: {exc}")

    run_tasks(
        tasks=tasks,
        worker=_load_spectra_from_mzml,
        workers=workers if len(file_tasks) > 1 else 1,
        on_start=_on_start,
        on_result=_on_result,
        on_error=_on_error,
    )

    return all_results


def load_spectra_streaming(
    variant_groups: Dict[str, VariantGroup],
    mzml_paths: List[Path],
    reference_substring: str = "",
    on_reference_loaded: Optional[Callable[[str, int], None]] = None,
    on_query_ready: Optional[Callable[[StreamingResult], None]] = None,
    on_file_start: Optional[Callable[[str, bool], None]] = None,
    on_file_complete: Optional[Callable[[str, int, int], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None,
    workers: int = 1,
    skip_ids: Optional[set] = None,
) -> Dict[str, List[LoadedSpectrum]]:
    """Load spectra with streaming: references first, then queries via callback.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup with matches.
        mzml_paths: List of available mzML paths.
        reference_substring: Substring (case-insensitive) identifying reference files.
        on_reference_loaded: Callback(variant_key, count) after references are loaded.
        on_query_ready: Callback(StreamingResult) for each query spectrum.
        on_file_start: Callback(filename, is_reference) when starting a file.
        on_file_complete: Callback(filename, loaded_count, not_found_count) when done.
        error_callback: Callback(error_message) for errors.
        workers: Number of parallel workers. 1 = sequential (default).
        skip_ids: Set of (variant_key, spectrum_id) tuples to skip (for resume).

    Returns:
        Dict of variant_key -> list of reference spectra.
    """
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}
    mzml_disambig = disambiguate_paths([str(p) for p in mzml_paths])

    reference_keys: Dict[str, Path] = {}
    query_keys: Dict[str, Path] = {}
    for path_key, path in mzml_path_map.items():
        if reference_substring and reference_substring.lower() in str(path).lower():
            reference_keys[path_key] = path
        else:
            query_keys[path_key] = path

    mzml_to_matches = build_mzml_to_matches_map(variant_groups, mzml_path_map)

    # Filter already-loaded spectra when resuming from partial cache
    if skip_ids:
        before = sum(len(m) for m in mzml_to_matches.values())
        mzml_to_matches = {
            pk: [(m, vk) for m, vk in ms if (vk, m.spectrum_id) not in skip_ids]
            for pk, ms in mzml_to_matches.items()
        }
        mzml_to_matches = {k: v for k, v in mzml_to_matches.items() if v}
        skipped = before - sum(len(m) for m in mzml_to_matches.values())
        if skipped > 0:
            logger.info(f"Resume: skipping {skipped} cached spectra, "
                        f"{sum(len(v) for v in mzml_to_matches.values())} remaining")

    reference_spectra: Dict[str, List[LoadedSpectrum]] = {k: [] for k in variant_groups}

    # Count queries per variant (for progress tracking)
    query_counts: Dict[str, int] = {k: 0 for k in variant_groups}
    for path_key, matches in mzml_to_matches.items():
        if path_key not in reference_keys:
            for _, variant_key in matches:
                query_counts[variant_key] += 1

    total_files = sum(1 for pk in list(reference_keys) + list(query_keys) if pk in mzml_to_matches)
    effective_workers = max(1, min(workers, total_files))

    # Phase 1: References
    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map, mzml_disambig)
    ref_results = _execute_load_tasks(ref_tasks, True, effective_workers,
                                      on_file_start, on_file_complete, error_callback)
    for res in ref_results:
        for variant_key, spec in res.loaded:
            reference_spectra[variant_key].append(spec)

    if on_reference_loaded:
        for variant_key, refs in reference_spectra.items():
            if refs:
                on_reference_loaded(variant_key, len(refs))

    logger.info(f"References loaded: {sum(len(v) for v in reference_spectra.values())} spectra")

    # Phase 2: Queries - stream via callback
    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map, mzml_disambig)
    query_results = _execute_load_tasks(query_tasks, False, effective_workers,
                                        on_file_start, on_file_complete, error_callback)

    query_indices: Dict[str, int] = {k: 0 for k in variant_groups}
    for res in query_results:
        for variant_key, loaded in res.loaded:
            if on_query_ready:
                variant_group = variant_groups[variant_key]
                on_query_ready(StreamingResult(
                    query_spectrum=loaded,
                    reference_spectra=reference_spectra[variant_key],
                    variant_key=variant_key,
                    base_sequence=variant_group.base_sequence,
                    modification_key=variant_group.modification_key,
                    query_index=query_indices[variant_key],
                    total_queries=query_counts[variant_key],
                ))
            query_indices[variant_key] += 1

    logger.info("Query streaming complete")
    return reference_spectra


def spectrum_lookup_key(spec: "LoadedSpectrum") -> str:
    """Return ``spectrum_id|peptidoform`` — stable PSM-unique key for lookup.

    The same physical spectrum may appear with different peptidoforms across
    PSMs, and spectrum IDs alone can collide across mzML files. Pair the key
    with ``build_psm_table()``'s ``psm_lookup`` dict to resolve ``psm_id`` when
    building peaks and sequence DataFrames.
    """
    return f"{spec.spectrum_id}|{spec.peptidoform}"
