"""Streaming mzML loader — orchestration + dispatch.

Coordinates reference-first ordering so downstream consumers can compute
query/reference metrics as each query arrives. The per-file worker itself
lives in ``core.spectra.mzml_workers``; parallel/sequential dispatch
lives in ``core.common.task_dispatch``.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from psm_utils import PSM, PSMList

from core.common.task_dispatch import run_tasks
from core.spectra.loaded_spectrum import LoadResult, LoadedSpectrum, StreamingResult
from core.spectra.mzml_workers import _load_spectra_from_mzml

logger = logging.getLogger(__name__)


def build_mzml_to_matches_map(
    variant_groups: Dict[str, PSMList],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Map each mzML path key to the list of ``(psm, variant_key)`` pairs to load.

    Prefers ``psm.metadata["mzml_resolved_path"]`` when set, then falls back to
    a stem-based lookup against ``psm.run`` (the canonical run name psm_utils
    parses from the idXML's ``spectra_data``).
    """
    # Build stem index: cheap reverse-lookup keyed on the mzML filename stem
    stem_to_key = {path.stem.lower(): key for key, path in available_mzml_paths.items()}

    # Resolver: prefer the PSM's pre-resolved path, fall back to stem match
    def resolve(psm: PSM) -> Optional[str]:
        resolved = psm.metadata.get("mzml_resolved_path", "")
        if resolved:
            key = str(Path(resolved).resolve()).lower()
            if key in available_mzml_paths:
                return key
        return stem_to_key.get(psm.run.lower())

    # Bucket: group each PSM under its resolved mzML path key
    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]] = {}
    for variant_key, psm_list in variant_groups.items():
        for psm in psm_list:
            if psm.run and (path_key := resolve(psm)):
                mzml_to_matches.setdefault(path_key, []).append((psm, variant_key))
    return mzml_to_matches


def _build_file_tasks(file_keys, mzml_to_matches, mzml_path_map):
    """Build ``(path_str, display_name, matches)`` tuples for files that have matches."""
    return [
        (str(mzml_path_map[key]), mzml_path_map[key].name, mzml_to_matches[key])
        for key in file_keys if key in mzml_to_matches
    ]


def _execute_load_tasks(file_tasks, is_reference, workers,
                        on_file_start, on_file_complete, error_callback) -> List[LoadResult]:
    """Load spectra from mzML files, sequentially or in parallel."""
    # Setup: short-circuit on empty, wire optional callbacks
    if not file_tasks:
        return []
    kind = "reference" if is_reference else "query"
    on_start = on_file_start or (lambda *_: None)
    on_complete = on_file_complete or (lambda *_: None)
    on_error = error_callback or (lambda *_: None)

    # Build & announce: pack is_reference into each task; fire start callbacks
    tasks = [(p, d, m, is_reference) for p, d, m in file_tasks]
    for _path, display_name, _matches, _is_ref in tasks:
        on_start(display_name, is_reference)

    # Dispatch & collect: surface worker errors, accumulate per-file results
    all_results: List[LoadResult] = []
    for task, ok, result in run_tasks(
        tasks=tasks, worker=_load_spectra_from_mzml, workers=workers,
    ):
        if not ok:
            on_error(f"Worker crashed loading {kind} {task[1]}: {result}")
            continue
        for error in result.errors:
            on_error(error)
        all_results.append(result)
        on_complete(result.display_name, result.loaded_count, result.not_found_count)
    return all_results


def _split_reference_and_query_keys(
    mzml_path_map: Dict[str, Path], reference_substring: str,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """Partition path-keys into ``(references, queries)`` by case-insensitive substring."""
    needle = reference_substring.lower() if reference_substring else None
    references: Dict[str, Path] = {}
    queries: Dict[str, Path] = {}
    for path_key, path in mzml_path_map.items():
        bucket = references if needle and needle in str(path).lower() else queries
        bucket[path_key] = path
    return references, queries


def _filter_already_loaded(
    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]],
    skip_ids: set,
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Drop any (variant_key, spectrum_id) pairs already present in *skip_ids*."""
    before = sum(len(matches) for matches in mzml_to_matches.values())
    filtered = {
        path_key: kept for path_key, matches in mzml_to_matches.items()
        if (kept := [(psm, vk) for psm, vk in matches if (vk, psm.spectrum_id) not in skip_ids])
    }
    after = sum(len(matches) for matches in filtered.values())
    if before > after:
        logger.info(f"Resume: skipping {before - after} cached spectra, {after} remaining")
    return filtered


def _collect_reference_spectra(
    ref_results: List[LoadResult], variant_groups: Dict[str, PSMList],
    on_reference_loaded: Optional[Callable[[str, int], None]],
) -> Dict[str, List[LoadedSpectrum]]:
    """Bucket reference spectra by variant_key and fire the per-variant callback."""
    # Bucket: initialise one slot per known variant, fan out the loaded spectra
    reference_spectra: Dict[str, List[LoadedSpectrum]] = {variant_key: [] for variant_key in variant_groups}
    for result in ref_results:
        for variant_key, spectrum in result.loaded:
            reference_spectra[variant_key].append(spectrum)

    # Notify: per-variant callback fires only when references were found
    notify = on_reference_loaded or (lambda *_: None)
    for variant_key, references in reference_spectra.items():
        if references:
            notify(variant_key, len(references))

    # Final summary
    logger.info(f"References loaded: {sum(len(r) for r in reference_spectra.values())} spectra")
    return reference_spectra


def _stream_query_results(
    query_results: List[LoadResult],
    reference_spectra: Dict[str, List[LoadedSpectrum]],
    on_query_ready: Optional[Callable[[StreamingResult], None]],
) -> None:
    """Fire ``on_query_ready`` once per loaded query spectrum."""
    if not on_query_ready:
        return
    for result in query_results:
        for variant_key, loaded in result.loaded:
            on_query_ready(StreamingResult(
                query_spectrum=loaded,
                reference_spectra=reference_spectra[variant_key],
                variant_key=variant_key,
            ))


def load_spectra_streaming(
    variant_groups: Dict[str, PSMList],
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
    """Stream-load references first, then queries; yields each query via ``on_query_ready``.

    ``reference_substring`` (case-insensitive) splits ``mzml_paths`` into
    references and queries. Callbacks fire on file start/complete/error plus
    per reference and per query. Returns ``{variant_key: [reference_spectra]}``.
    """
    # Partition inputs: split paths into reference/query buckets, build PSM map
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}
    reference_keys, query_keys = _split_reference_and_query_keys(mzml_path_map, reference_substring)

    mzml_to_matches = build_mzml_to_matches_map(variant_groups, mzml_path_map)
    if skip_ids:
        mzml_to_matches = _filter_already_loaded(mzml_to_matches, skip_ids)

    # References first: load + bucket so queries can compute metrics immediately
    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map)
    ref_results = _execute_load_tasks(
        ref_tasks, True, workers, on_file_start, on_file_complete, error_callback,
    )
    reference_spectra = _collect_reference_spectra(ref_results, variant_groups, on_reference_loaded)

    # Queries next: stream each loaded query through on_query_ready
    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map)
    query_results = _execute_load_tasks(
        query_tasks, False, workers, on_file_start, on_file_complete, error_callback,
    )
    _stream_query_results(query_results, reference_spectra, on_query_ready)

    logger.info("Query streaming complete")
    return reference_spectra
