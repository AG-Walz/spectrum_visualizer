"""Streaming mzML loader — orchestration + dispatch.

Coordinates reference-first ordering so downstream consumers can compute
query/reference metrics as each query arrives. The per-file worker itself
lives in ``mspecviz.spectra.loader_workers``; parallel/sequential dispatch
lives in ``mspecviz.common.parallel``.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from psm_utils import PSM

from mspecviz.common.parallel import run_tasks
from mspecviz.search.models import VariantGroup
from mspecviz.spectra.loader_models import LoadResult, LoadedSpectrum, StreamingResult
from mspecviz.spectra.loader_workers import _load_spectra_from_mzml

logger = logging.getLogger(__name__)


def build_mzml_to_matches_map(
    variant_groups: Dict[str, VariantGroup],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Map each mzML path key to the list of ``(psm, variant_key)`` pairs to load.

    Prefers ``psm.metadata["mzml_resolved_path"]`` when set, then falls back to
    a stem-based lookup against ``psm.run`` (the canonical run name psm_utils
    parses from the idXML's ``spectra_data``). The unique-basename invariant
    enforced at input makes the fallback unambiguous.
    """
    stem_to_key: Dict[str, str] = {
        path.stem.lower(): key for key, path in available_mzml_paths.items()
    }

    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]] = {}

    for variant_key, variant_group in variant_groups.items():
        for psm in variant_group.matches:
            if not psm.run:
                continue

            path_key = None
            resolved = psm.metadata.get("mzml_resolved_path", "")
            if resolved:
                resolved_key = str(Path(resolved).resolve()).lower()
                if resolved_key in available_mzml_paths:
                    path_key = resolved_key

            if path_key is None:
                path_key = stem_to_key.get(psm.run.lower())

            if path_key:
                mzml_to_matches.setdefault(path_key, []).append((psm, variant_key))

    return mzml_to_matches


def _build_file_tasks(file_keys, mzml_to_matches, mzml_path_map):
    """Build (path_str, display_name, matches) tuples for files that have matches."""
    tasks = []
    for path_key in file_keys:
        if path_key not in mzml_to_matches:
            continue
        mzml_path = mzml_path_map[path_key]
        tasks.append((str(mzml_path), mzml_path.name, mzml_to_matches[path_key]))
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

    # Mirror the original eager-callback behaviour: announce all files up front.
    if on_file_start:
        for _path, display_name, _matches, _is_ref in tasks:
            on_file_start(display_name, is_reference)

    for task, ok, result in run_tasks(
        tasks=tasks,
        worker=_load_spectra_from_mzml,
        workers=workers if len(file_tasks) > 1 else 1,
    ):
        if not ok:
            display_name = task[1]
            if error_callback:
                error_callback(f"Worker crashed loading {kind} {display_name}: {result}")
            continue
        if error_callback:
            for error in result.errors:
                error_callback(error)
        all_results.append(result)
        if on_file_complete:
            on_file_complete(result.display_name, result.loaded_count, result.not_found_count)

    return all_results


def _split_reference_and_query_keys(
    mzml_path_map: Dict[str, Path], reference_substring: str,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """Partition path-keys into ``(references, queries)`` by case-insensitive substring."""
    reference_keys: Dict[str, Path] = {}
    query_keys: Dict[str, Path] = {}
    needle = reference_substring.lower() if reference_substring else None
    for path_key, path in mzml_path_map.items():
        if needle and needle in str(path).lower():
            reference_keys[path_key] = path
        else:
            query_keys[path_key] = path
    return reference_keys, query_keys


def _filter_already_loaded(
    mzml_to_matches: Dict[str, List[Tuple[PSM, str]]],
    skip_ids: set,
) -> Dict[str, List[Tuple[PSM, str]]]:
    """Drop any (variant_key, spectrum_id) pairs already present in *skip_ids*."""
    before = sum(len(matches) for matches in mzml_to_matches.values())
    filtered = {
        path_key: [(psm, variant_key) for psm, variant_key in matches
                   if (variant_key, psm.spectrum_id) not in skip_ids]
        for path_key, matches in mzml_to_matches.items()
    }
    filtered = {path_key: matches for path_key, matches in filtered.items() if matches}
    skipped = before - sum(len(matches) for matches in filtered.values())
    if skipped > 0:
        logger.info(
            f"Resume: skipping {skipped} cached spectra, "
            f"{sum(len(matches) for matches in filtered.values())} remaining"
        )
    return filtered


def _collect_reference_spectra(
    ref_results: List[LoadResult], variant_groups: Dict[str, VariantGroup],
    on_reference_loaded: Optional[Callable[[str, int], None]],
) -> Dict[str, List[LoadedSpectrum]]:
    """Bucket reference spectra by variant_key and fire the per-variant callback."""
    reference_spectra: Dict[str, List[LoadedSpectrum]] = {variant_key: [] for variant_key in variant_groups}
    for result in ref_results:
        for variant_key, spectrum in result.loaded:
            reference_spectra[variant_key].append(spectrum)
    if on_reference_loaded:
        for variant_key, references in reference_spectra.items():
            if references:
                on_reference_loaded(variant_key, len(references))
    logger.info(
        f"References loaded: "
        f"{sum(len(refs) for refs in reference_spectra.values())} spectra"
    )
    return reference_spectra


def _stream_query_results(
    query_results: List[LoadResult],
    reference_spectra: Dict[str, List[LoadedSpectrum]],
    variant_groups: Dict[str, VariantGroup],
    on_query_ready: Optional[Callable[[StreamingResult], None]],
) -> None:
    """Fire ``on_query_ready`` once per loaded query spectrum."""
    if not on_query_ready:
        return
    for result in query_results:
        for variant_key, loaded in result.loaded:
            variant_group = variant_groups[variant_key]
            on_query_ready(StreamingResult(
                query_spectrum=loaded,
                reference_spectra=reference_spectra[variant_key],
                variant_key=variant_key,
                base_sequence=variant_group.base_sequence,
                modification_key=variant_group.modification_key,
            ))


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
    """Stream-load references first, then queries; yields each query via ``on_query_ready``.

    ``reference_substring`` (case-insensitive) splits ``mzml_paths`` into
    references and queries. Callbacks fire on file start/complete/error plus
    per reference and per query. Returns ``{variant_key: [reference_spectra]}``.
    """
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}
    reference_keys, query_keys = _split_reference_and_query_keys(
        mzml_path_map, reference_substring,
    )

    mzml_to_matches = build_mzml_to_matches_map(variant_groups, mzml_path_map)
    if skip_ids:
        mzml_to_matches = _filter_already_loaded(mzml_to_matches, skip_ids)

    total_files = sum(
        1 for pk in list(reference_keys) + list(query_keys) if pk in mzml_to_matches
    )
    effective_workers = max(1, min(workers, total_files))

    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map)
    ref_results = _execute_load_tasks(
        ref_tasks, True, effective_workers,
        on_file_start, on_file_complete, error_callback,
    )
    reference_spectra = _collect_reference_spectra(
        ref_results, variant_groups, on_reference_loaded,
    )

    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map)
    query_results = _execute_load_tasks(
        query_tasks, False, effective_workers,
        on_file_start, on_file_complete, error_callback,
    )
    _stream_query_results(query_results, reference_spectra, variant_groups, on_query_ready)

    logger.info("Query streaming complete")
    return reference_spectra
