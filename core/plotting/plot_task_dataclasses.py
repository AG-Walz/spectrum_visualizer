"""Shared, picklable dataclasses that carry one mirror/single plot job — and
its result — across the process boundary into the parallel render workers.
ComparisonTask is the input work unit; ComparisonResult is the output.
Both are deliberately plain dataclasses holding only picklable values so they
can be shipped to a ProcessPoolExecutor worker (built by
core.pipeline.query_processing.build_comparison_tasks, consumed by
core.plotting.plot_worker.process_single_task; re-exported via
core.plotting.__init__).
"""
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ComparisonTask:
    """Picklable bundle of everything a worker needs to render one plot.
    One instance describes a single mirror or single-spectrum plot job:
    the annotated spectra, the labels, the on-disk output locations, and the
    export parameters. Created in batches by
    core.pipeline.query_processing.build_comparison_tasks and consumed by
    core.plotting.plot_worker.process_single_task inside a worker process.
    All fields hold picklable values only (DataFrames, strings, numbers) so the
    whole object survives being shipped to a ProcessPoolExecutor worker.
    """
    task_id: int          # Stable index — paired back to a ComparisonResult
    task_type: str        # "mirror" (query vs reference) or "single" (query only)

    # --- Annotated spectra -------------------------------------------------
    # Annotated peak DataFrames carried by value so workers need no shared state.
    query_frame: pd.DataFrame
    reference_frame: Optional[pd.DataFrame] = None  # None for "single" tasks

    # --- Plot labels -------------------------------------------------------
    query_label: str = ""
    reference_label: str = ""
    plot_title: str = ""

    # --- Export paths ------------------------------------------------------
    # Plain strings — workers create any missing subdirectories themselves.
    folder_name: str = ""
    safe_filename: str = ""
    png_dir: str = ""
    svg_dir: str = ""
    pdf_dir: str = ""

    # --- Plot / export parameters -----------------------------------------
    export_formats: list[str] = field(default_factory=lambda: ["png"])
    figure_size: tuple[int, int] = (12, 8)  # inches (width, height)
    export_resolution: int = 300            # DPI for raster output

    # Output path for the interactive Plotly HTML — empty string skips HTML.
    html_plot_path: str = ""

    # Skip plot rendering entirely (metrics-only run — plot already on disk).
    skip_plot: bool = False

    # Pre-formatted metrics box: lets the caller compute metrics once in the
    # main process and reuse the text across workers instead of recomputing.
    metrics_text: str = ""

    # --- Summary metadata --------------------------------------------------
    # Echoed verbatim into the summary-CSV row built for mirror tasks.
    query_file: str = ""
    reference_file: str = ""
    peptide: str = ""
    modification: str = ""
    charge: int = 0


@dataclass
class ComparisonResult:
    """Worker output for one ComparisonTask.
    Returned by core.plotting.plot_worker.process_single_task and collected
    back in the main process to report progress and assemble the summary CSV.
    Carries the success flag, the on-disk path written for each rendered format,
    and — for mirror tasks — the metrics row destined for the summary CSV.
    """
    task_id: int        # Matches the originating ComparisonTask.task_id
    task_type: str       # Copied from the task ("mirror" or "single")
    success: bool        # False when rendering raised — see error
    error: str = ""      # Exception text when success is False, else ""

    # --- Rendered output paths --------------------------------------------
    # Empty string for any format not requested or not produced.
    png_path: str = ""
    svg_path: str = ""
    pdf_path: str = ""
    html_path: str = ""

    # Per-task metrics row for the summary CSV — populated for mirror tasks only.
    summary_row: Optional[dict] = None
