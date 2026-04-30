"""Streaming mzML loader — orchestration + dispatch.

Coordinates reference-first ordering so downstream consumers can compute
query/reference metrics as each query arrives. The per-file worker itself
lives in ``mspecviz.spectra.loader_workers``; parallel/sequential dispatch
lives in ``mspecviz.common.parallel``.
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from mspecviz.common.parallel import run_tasks
from mspecviz.search.models import PeptideMatch, VariantGroup
from mspecviz.spectra.loader_models import LoadResult, LoadedSpectrum, StreamingResult
from mspecviz.spectra.loader_workers import _load_spectra_from_mzml

logger = logging.getLogger(__name__)


def build_mzml_to_matches_map(
    variant_groups: Dict[str, VariantGroup],
    available_mzml_paths: Dict[str, Path],
) -> Dict[str, List[Tuple[PeptideMatch, str]]]:
    """Map each mzML path key to the list of (match, variant_key) pairs to load.

    Prefers ``mzml_resolved_path`` when set, then falls back to a bare-name
    lookup. The unique-basename invariant (enforced at input) guarantees the
    fallback is unambiguous.
    """
    name_to_key: Dict[str, str] = {
        path.name.lower(): key for key, path in available_mzml_paths.items()
    }

    mzml_to_matches: Dict[str, List[Tuple[PeptideMatch, str]]] = {}

    for variant_key, variant_group in variant_groups.items():
        for match in variant_group.matches:
            if not match.mzml_filename:
                continue

            path_key = None
            if match.mzml_resolved_path:
                resolved_key = str(Path(match.mzml_resolved_path).resolve()).lower()
                if resolved_key in available_mzml_paths:
                    path_key = resolved_key

            if path_key is None:
                path_key = name_to_key.get(match.mzml_filename.lower())

            if path_key:
                mzml_to_matches.setdefault(path_key, []).append((match, variant_key))

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
    """Stream-load references first, then queries; yields each query via ``on_query_ready``.

    ``reference_substring`` (case-insensitive) splits ``mzml_paths`` into
    references and queries. Callbacks fire on file start/complete/error plus
    per reference and per query. Returns ``{variant_key: [reference_spectra]}``.
    """
    mzml_path_map = {str(p.resolve()).lower(): p for p in mzml_paths}

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

    total_files = sum(1 for pk in list(reference_keys) + list(query_keys) if pk in mzml_to_matches)
    effective_workers = max(1, min(workers, total_files))

    # Phase 1: References
    ref_tasks = _build_file_tasks(reference_keys, mzml_to_matches, mzml_path_map)
    ref_results = _execute_load_tasks(ref_tasks, True, effective_workers,
                                      on_file_start, on_file_complete, error_callback)
    for res in ref_results:
        for variant_key, spec in res.loaded:
            reference_spectra[variant_key].append(spec)

    if on_reference_loaded:
        for variant_key, references in reference_spectra.items():
            if references:
                on_reference_loaded(variant_key, len(references))

    logger.info(f"References loaded: {sum(len(v) for v in reference_spectra.values())} spectra")

    # Phase 2: Queries - stream via callback
    query_tasks = _build_file_tasks(query_keys, mzml_to_matches, mzml_path_map)
    query_results = _execute_load_tasks(query_tasks, False, effective_workers,
                                        on_file_start, on_file_complete, error_callback)

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
                ))

    logger.info("Query streaming complete")
    return reference_spectra
