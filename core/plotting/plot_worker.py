"""ProcessPoolExecutor-picklable worker that renders one peptide-spectrum comparison to disk.
Role — leaf of the export pipeline: a stateless, picklable worker function executed
inside a ProcessPoolExecutor pool, one invocation per ComparisonTask.
Produces — for a single task, computes the 7 comparison metrics (mirror tasks only)
and writes the requested output files: a matplotlib static figure as PNG/SVG/PDF
under per-variant subdirectories, and optionally a standalone Plotly HTML plot.
Each invocation returns a ComparisonResult carrying the written paths plus one
summary-CSV row, never raising — failures are reported via the result's error field.
Call relationships — process_single_task is dispatched by core.export.export_pipeline
through core.helpers.run_tasks across the worker pool (re-exported via
core.plotting.process_single_task). It calls core.annotation.comparison_metrics
for metric computation and formatting, and core.plotting.mirror.create_spectrum_plot
for both the matplotlib and Plotly figures. ComparisonTask / ComparisonResult
come from core.plotting.plot_task_dataclasses.
"""
import logging
from pathlib import Path

# Force the non-interactive Agg backend before pyplot is imported — worker processes
# have no display and must render figures purely to disk.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from core.annotation.comparison_metrics import compute_all_metrics, format_metrics_text
from core.plotting.mirror import create_spectrum_plot
from core.plotting.plot_task_dataclasses import ComparisonResult, ComparisonTask

logger = logging.getLogger(__name__)


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
    Runs inside a ProcessPoolExecutor worker: the argument and the returned
    ComparisonResult must stay picklable, so only plain data crosses the
    process boundary. Computes the comparison metrics (mirror tasks only),
    renders the static matplotlib figure into each requested format, and
    optionally writes a Plotly HTML plot.
    Returns a ComparisonResult with success=True and the written file
    paths on completion. Any exception is caught and converted into a result
    with success=False and the message in error — the worker never
    raises, so one bad task cannot abort the whole pool.
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
    Each format lands in its own per-variant subdirectory
    (<format_dir>/<folder_name>/<safe_filename>.<ext>), created on demand.
    Returns the triple (png_path, svg_path, pdf_path) — a format is skipped
    (leaving an empty string) when it is absent from task.export_formats or
    its destination directory is unset.
    """
    # Pair each format with its configured destination directory.
    format_dirs = {"png": task.png_dir, "svg": task.svg_dir, "pdf": task.pdf_dir}
    paths = {"png": "", "svg": "", "pdf": ""}

    for output_format, format_dir in format_dirs.items():
        # Skip formats the task did not request or that have no target directory.
        if output_format not in task.export_formats or not format_dir:
            continue
        # Variants get their own subfolder so plots stay grouped per peptide/reference.
        variant_dir = Path(format_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        output_path = variant_dir / f"{task.safe_filename}.{output_format}"
        save_kwargs: dict = {"bbox_inches": "tight", "facecolor": "white", "format": output_format}
        # DPI only meaningfully applies to the raster PNG; SVG/PDF are vector.
        if output_format == "png":
            save_kwargs["dpi"] = task.export_resolution
        fig.savefig(str(output_path), **save_kwargs)
        paths[output_format] = str(output_path)

    return paths["png"], paths["svg"], paths["pdf"]


def _save_task_html(task: ComparisonTask) -> str:
    """Render the Plotly mirror/single plot and write it as standalone HTML.
    Returns the written path on success, or "" if an OSError occurred
    while creating the directory or writing the file — the failure is logged
    as a warning rather than propagated, so a single unwritable HTML target
    does not fail the whole task.
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
