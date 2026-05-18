"""Top-level streaming entry: load spectra, annotate each PSM, buffer results (called by ``cli.main`` and ``streamlit_app.views.search_load``)."""
from pathlib import Path
from typing import Callable, Optional

from core.pipeline.query_processing import process_query_result
from core.spectra import load_spectra_streaming


def make_query_ready_callback(
    *,
    psm_lookup: dict,
    annotation_kwargs: dict,
    reference_store: dict,
    koina_kwargs: dict,
    annotated_queries: list,
    comparison_rows: list,
    on_progress: Callable[[int], None],
):
    """Build the per-PSM callback that annotates, appends to buffers, and ticks ``on_progress``."""
    loaded_count = [0]

    def on_query_ready(result):
        processed = process_query_result(
            result, psm_lookup, annotation_kwargs, reference_store, **koina_kwargs,
        )
        if processed is None:
            return
        psm_id, query_spectrum, query_frame, rows = processed
        annotated_queries.append((psm_id, query_spectrum, query_frame))
        comparison_rows.extend(rows)
        loaded_count[0] += 1
        on_progress(loaded_count[0])

    return on_query_ready


def run_streaming_annotation(
    *,
    variant_groups: dict,
    mzml_paths,
    reference_substring: str,
    psm_lookup: dict,
    annotation_kwargs: dict,
    koina_kwargs: dict,
    workers: int,
    on_progress: Optional[Callable[[int], None]] = None,
    on_file_start: Optional[Callable] = None,
    on_file_complete: Optional[Callable] = None,
    on_error: Optional[Callable[[str], None]] = None,
) -> tuple[list, dict, list]:
    """Stream-load and annotate every PSM; return ``(annotated_queries, reference_store, comparison_rows)``."""
    # Init buffers: shared accumulators populated by the on_query_ready callback
    annotated_queries: list = []
    reference_store: dict = {}
    comparison_rows: list = []
    noop = lambda *_: None

    # Build callback: per-PSM closure that annotates and appends to buffers
    on_query_ready = make_query_ready_callback(
        psm_lookup=psm_lookup,
        annotation_kwargs=annotation_kwargs,
        reference_store=reference_store,
        koina_kwargs=koina_kwargs,
        annotated_queries=annotated_queries,
        comparison_rows=comparison_rows,
        on_progress=on_progress or noop,
    )

    # Dispatch: stream every mzML, fanning lifecycle callbacks to the frontend
    load_spectra_streaming(
        variant_groups=variant_groups,
        mzml_paths=[Path(path) for path in mzml_paths],
        reference_substring=reference_substring,
        on_reference_loaded=noop,
        on_query_ready=on_query_ready,
        on_file_start=on_file_start or noop,
        on_file_complete=on_file_complete or noop,
        error_callback=on_error or noop,
        workers=workers,
    )

    return annotated_queries, reference_store, comparison_rows
