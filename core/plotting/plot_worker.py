"""ProcessPoolExecutor-picklable worker that renders one peptide-spectrum comparison to disk.

This module bundles the picklable work-unit dataclasses (`ComparisonTask` input,
`ComparisonResult` output) and the stateless worker entry point
(`process_single_task`) that runs inside a ProcessPoolExecutor pool.

For a single task, the worker computes the 7 comparison metrics (mirror tasks
only) and writes the requested output files: a Matplotlib static figure as
PNG/SVG/PDF under per-peptidoform subdirectories, and optionally a standalone
Plotly HTML plot. Each invocation returns a ComparisonResult carrying the
written paths plus one summary-CSV row, never raising — failures are reported
via the result's error field so one bad task cannot abort the whole pool.

`process_single_task` is dispatched by core.export.run_export through
core.helpers.run_tasks across the worker pool (re-exported via
core.plotting.process_single_task). It calls core.comparison_metrics
for metric computation and formatting, and core.plotting.mirror.create_spectrum_plot
for both the Matplotlib and Plotly figures.
"""
import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Force the non-interactive Agg backend before pyplot is imported — worker processes
# have no display and must render figures purely to disk.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from core.comparison_metrics import compute_all_metrics, format_metrics_text
from core.plotting.mirror import create_spectrum_plot

logger = logging.getLogger(__name__)


# ============================================================================
# Picklable work-unit dataclasses (shipped across the ProcessPool boundary)
# ============================================================================

@dataclass
class ComparisonTask:
    """Picklable bundle of everything a worker needs to render one plot.
    Describes a single mirror or single-spectrum plot job: the annotated spectra, the labels,
    the on-disk output locations, and the export parameters. Built in batches by build_plot_tasks
    and consumed by process_single_task inside a worker process; all fields hold picklable values
    only (DataFrames, strings, numbers) so the whole object survives the ProcessPoolExecutor boundary.
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
    Returned by process_single_task and collected in the main process to report
    progress and assemble the summary CSV. Carries the success flag, the path
    written per format, and — for mirror tasks — the summary-CSV metrics row.
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


# ============================================================================
# Worker entry point + per-format save helpers
# ============================================================================

def _build_summary_row(task: ComparisonTask, metrics: dict) -> dict:
    """Compose one row of the comparison CSV for a single mirror task.
    Returns a flat dict keying the identity columns (query_file,
    ref_file, peptide, modification, charge) followed by every
    computed metric spread in from metrics — one record per output row.
    """
    return {
        "query_file": task.query_file,
        "ref_file": task.reference_file,
        "peptide": task.peptide,
        "modification": task.modification,
        "charge": task.charge,
        **metrics,
    }


def process_single_task(task: ComparisonTask) -> ComparisonResult:
    """Worker entry point — render one ComparisonTask to disk.
    Runs inside a ProcessPoolExecutor worker, so the argument and the returned ComparisonResult
    must stay picklable. Computes the comparison metrics (mirror tasks only), renders the static
    matplotlib figure into each requested format, and optionally writes a Plotly HTML plot. Returns
    a ComparisonResult with success=True on completion; any exception becomes success=False, never raised.
    """
    try:
        # Classify: mirror tasks need a reference; remember slots for output paths
        is_mirror = task.task_type == "mirror" and task.reference_frame is not None
        summary_row = None
        png_path = svg_path = pdf_path = html_path = ""

        # Metrics: always computed for mirror tasks, even when the plot is skipped
        metrics = None
        if is_mirror:
            metrics = compute_all_metrics(task.query_frame, task.reference_frame)
            summary_row = _build_summary_row(task, metrics)

        # Static figure: matplotlib render, then write every requested format
        if not task.skip_plot:
            fig = create_spectrum_plot(
                task.query_frame,
                task.reference_frame if is_mirror else None,
                backend="ms_matplotlib",
                title=task.plot_title,
                query_label=task.query_label,
                reference_label=task.reference_label,
                metrics_text=(task.metrics_text or format_metrics_text(metrics)) if is_mirror else "",
                figure_size=task.figure_size,
            )
            png_path, svg_path, pdf_path = _save_task_figure(fig, task)
            plt.close(fig)

        # Interactive HTML: only when the task requests a Plotly export
        if task.html_plot_path:
            html_path = _save_task_html(task)

        # Return: hand the populated result back to the orchestrator
        return ComparisonResult(
            task_id=task.task_id, task_type=task.task_type, success=True,
            png_path=png_path, svg_path=svg_path, pdf_path=pdf_path,
            html_path=html_path, summary_row=summary_row,
        )
    except Exception as error:
        return ComparisonResult(
            task_id=task.task_id, task_type=task.task_type,
            success=False, error=str(error),
        )


def _save_task_figure(fig, task: ComparisonTask) -> tuple[str, str, str]:
    """Write the rendered matplotlib figure to each requested file format.
    Each format lands in its own per-peptidoform subdirectory
    (<format_dir>/<folder_name>/<safe_filename>.<ext>), created on demand. Returns
    (png_path, svg_path, pdf_path); a skipped format leaves an empty string.
    """
    # Pair each format with its configured destination directory.
    format_dirs = {"png": task.png_dir, "svg": task.svg_dir, "pdf": task.pdf_dir}
    paths = {"png": "", "svg": "", "pdf": ""}

    for output_format, format_dir in format_dirs.items():
        # Skip formats the task did not request or that have no target directory.
        if output_format not in task.export_formats or not format_dir:
            continue
        # Peptidoforms get their own subfolder so plots stay grouped per peptide/reference.
        peptidoform_dir = Path(format_dir) / task.folder_name
        peptidoform_dir.mkdir(parents=True, exist_ok=True)
        output_path = peptidoform_dir / f"{task.safe_filename}.{output_format}"
        save_kwargs: dict = {"bbox_inches": "tight", "facecolor": "white", "format": output_format}
        # DPI only meaningfully applies to the raster PNG; SVG/PDF are vector.
        if output_format == "png":
            save_kwargs["dpi"] = task.export_resolution
        # Atomic write: render to a unique temp file in the same directory, then os.replace
        # into place. Two parallel workers that resolve to the same output_path can otherwise
        # interleave their bytes and leave a valid-header/garbled-body file (corrupt PNG IDAT);
        # os.replace is atomic within one filesystem, so every file on disk stays complete.
        tmp_fd, tmp_name = tempfile.mkstemp(dir=peptidoform_dir, suffix=f".{output_format}.tmp")
        os.close(tmp_fd)
        fig.savefig(tmp_name, **save_kwargs)
        os.replace(tmp_name, output_path)
        paths[output_format] = str(output_path)

    return paths["png"], paths["svg"], paths["pdf"]


def _save_task_html(task: ComparisonTask) -> str:
    """Render the Plotly mirror/single plot and write it as standalone HTML.
    Returns the written path, or "" on OSError (logged as a warning rather than
    propagated, so a single unwritable HTML target does not fail the whole task).
    """
    try:
        out = Path(task.html_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        # A mirror plot needs both a query and a reference frame; otherwise render single.
        is_mirror = task.task_type == "mirror" and task.reference_frame is not None
        fig = create_spectrum_plot(
            task.query_frame,
            task.reference_frame if is_mirror else None,
            backend="ms_plotly",
            title=task.plot_title,
            query_label=task.query_label,
            reference_label=task.reference_label,
            metrics_text=task.metrics_text if is_mirror else "",
        )
        # include_plotlyjs="cdn" keeps the file small — the browser pulls plotly.js from a CDN.
        fig.write_html(str(out), include_plotlyjs="cdn")
        return str(out)

    except OSError as error:
        logger.warning(f"HTML plot save failed for task {task.task_id}: {error}")
        return ""
