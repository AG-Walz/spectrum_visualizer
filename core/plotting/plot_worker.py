"""ProcessPoolExecutor-picklable worker: compute comparison metrics + render mirror/single plots to disk (dispatched by ``core.export.export_pipeline``)."""
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from core.annotation.comparison_metrics import compute_all_metrics, format_metrics_text
from core.plotting.mirror import create_spectrum_plot
from core.plotting.plot_task_dataclasses import ComparisonResult, ComparisonTask

logger = logging.getLogger(__name__)


def _build_summary_row(task: ComparisonTask, metrics: dict) -> dict:
    """Compose one row of the comparison CSV (file names, peptide identity, charge, plus every computed metric)."""
    return {
        "query_file": task.query_file,
        "ref_file": task.reference_file,
        "peptide": task.peptide,
        "modification": task.modification,
        "charge": task.charge,
        **metrics,
    }


def process_single_task(task: ComparisonTask) -> ComparisonResult:
    """Worker entry point: compute metrics for mirror tasks, render the figure, and write each requested output format to disk."""
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
    """Write the matplotlib figure to each requested format under per-variant subdirs; return ``(png_path, svg_path, pdf_path)`` with empty strings for skipped formats."""
    format_dirs = {"png": task.png_dir, "svg": task.svg_dir, "pdf": task.pdf_dir}
    paths = {"png": "", "svg": "", "pdf": ""}

    for output_format, format_dir in format_dirs.items():
        if output_format not in task.export_formats or not format_dir:
            continue
        variant_dir = Path(format_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        output_path = variant_dir / f"{task.safe_filename}.{output_format}"
        save_kwargs: dict = {"bbox_inches": "tight", "facecolor": "white", "format": output_format}
        if output_format == "png":
            save_kwargs["dpi"] = task.export_resolution
        fig.savefig(str(output_path), **save_kwargs)
        paths[output_format] = str(output_path)

    return paths["png"], paths["svg"], paths["pdf"]


def _save_task_html(task: ComparisonTask) -> str:
    """Render the Plotly mirror/single plot and write it as standalone HTML; return the output path, or ``""`` if the write failed."""
    try:
        out = Path(task.html_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)
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
        fig.write_html(str(out), include_plotlyjs="cdn")
        return str(out)

    except OSError as error:
        logger.warning(f"HTML plot save failed for task {task.task_id}: {error}")
        return ""
