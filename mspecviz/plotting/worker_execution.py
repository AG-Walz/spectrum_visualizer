"""Parallel plot generation via ProcessPoolExecutor.

Main process annotates spectra, worker pool runs ``process_single_task()``
(metrics + matplotlib/plotly), main process collects ``ComparisonResults``.
Task/result dataclasses live in ``worker_tasks`` so lightweight callers can
construct tasks without pulling matplotlib.
"""
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mspecviz.annotation.metrics_set import compute_all_metrics, format_metrics_text
from mspecviz.plotting.mirror_matplotlib import (
    create_mirror_plot,
    create_single_spectrum_plot,
)
from mspecviz.plotting.mirror_plotly import (
    create_plotly_mirror_plot,
    create_plotly_single_plot,
)
from mspecviz.plotting.worker_tasks import ComparisonResult, ComparisonTask

logger = logging.getLogger(__name__)


def _build_summary_row(task: ComparisonTask, metrics: dict) -> dict:
    """Compose the per-task row that becomes one line of the comparison CSV."""
    return {
        "query_file": task.query_file,
        "ref_file": task.reference_file,
        "peptide": task.peptide,
        "modification": task.modification,
        "charge": task.charge,
        **metrics,
    }


def _render_mirror_figure(task: ComparisonTask, metrics: dict):
    """Build the matplotlib mirror figure for a task (caller closes the figure)."""
    metrics_text = task.metrics_text or format_metrics_text(metrics)
    return create_mirror_plot(
        upper_spectrum=task.query_frame,
        lower_spectrum=task.reference_frame,
        upper_label=task.query_label,
        lower_label=task.reference_label,
        title=task.plot_title,
        metrics_text=metrics_text,
        figure_size=task.figure_size,
        export_resolution=task.export_resolution,
        normalize=task.normalize,
    )


def _render_single_figure(task: ComparisonTask):
    """Build the matplotlib single-spectrum figure for a task."""
    return create_single_spectrum_plot(
        spectrum=task.query_frame,
        title=task.plot_title,
        label=task.query_label,
        figure_size=task.figure_size,
        export_resolution=task.export_resolution,
        normalize=task.normalize,
    )


def _save_static_and_close(fig, task: ComparisonTask) -> Tuple[str, str, str]:
    """Save *fig* via ``_save_task_figure`` then close + delete it."""
    png_path, svg_path, pdf_path = _save_task_figure(fig, task)
    plt.close(fig)
    del fig
    return png_path, svg_path, pdf_path


def _process_mirror_task(task: ComparisonTask) -> Tuple[dict, Tuple[str, str, str], str]:
    """Run the mirror branch: metrics + (optional) static + (optional) HTML."""
    metrics = compute_all_metrics(task.query_frame, task.reference_frame)
    summary_row = _build_summary_row(task, metrics)
    paths = ("", "", "")
    if not task.skip_plot:
        fig = _render_mirror_figure(task, metrics)
        paths = _save_static_and_close(fig, task)
    html_path = _save_task_html(task) if task.html_plot_path else ""
    return summary_row, paths, html_path


def _process_single_task_branch(task: ComparisonTask) -> Tuple[Tuple[str, str, str], str]:
    """Run the single-spectrum branch: optional static + optional HTML."""
    paths = ("", "", "")
    if not task.skip_plot:
        fig = _render_single_figure(task)
        paths = _save_static_and_close(fig, task)
    html_path = _save_task_html(task) if task.html_plot_path else ""
    return paths, html_path


def process_single_task(task: ComparisonTask) -> ComparisonResult:
    """Worker: compute metrics, create figure, save to disk. Skips plot if skip_plot is set."""
    try:
        summary_row = None
        png_path = svg_path = pdf_path = html_path = ""

        if task.task_type == "mirror" and task.reference_frame is not None:
            summary_row, (png_path, svg_path, pdf_path), html_path = _process_mirror_task(task)
        elif not task.skip_plot:
            (png_path, svg_path, pdf_path), html_path = _process_single_task_branch(task)

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


# Format -> savefig keyword and extension mapping
_FORMAT_SAVE_KWARGS: Dict[str, dict] = {
    "png": {"format": "png"},
    "svg": {"format": "svg"},
    "pdf": {"format": "pdf"},
}


def _save_task_figure(fig, task: ComparisonTask) -> Tuple[str, str, str]:
    """Save a matplotlib figure to each requested format. Returns (png_path, svg_path, pdf_path)."""
    format_dirs = {"png": task.png_dir, "svg": task.svg_dir, "pdf": task.pdf_dir}
    paths = {"png": "", "svg": "", "pdf": ""}

    for output_format, format_dir in format_dirs.items():
        if output_format not in task.export_formats or not format_dir:
            continue
        variant_dir = Path(format_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        output_path = variant_dir / f"{task.safe_filename}.{output_format}"
        save_kwargs = {"bbox_inches": "tight", "facecolor": "white",
                       **_FORMAT_SAVE_KWARGS[output_format]}
        if output_format == "png":
            save_kwargs["dpi"] = task.export_resolution
        fig.savefig(str(output_path), **save_kwargs)
        paths[output_format] = str(output_path)

    return paths["png"], paths["svg"], paths["pdf"]


def _save_task_html(task: ComparisonTask) -> str:
    """Save interactive Plotly HTML plot. Returns path or "" on failure."""
    if not task.html_plot_path:
        return ""

    try:
        out = Path(task.html_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)

        if task.task_type == "mirror" and task.reference_frame is not None:
            fig = create_plotly_mirror_plot(
                upper_spectrum=task.query_frame,
                lower_spectrum=task.reference_frame,
                title=task.plot_title,
                normalize=task.normalize,
            )
        else:
            fig = create_plotly_single_plot(
                spectrum=task.query_frame,
                title=task.plot_title,
                normalize=task.normalize,
            )

        fig.write_html(str(out), include_plotlyjs="cdn")
        return str(out)

    except OSError as error:
        logger.warning(f"HTML plot save failed for task {task.task_id}: {error}")
        return ""


def run_parallel_comparisons(
    tasks: List[ComparisonTask],
    workers: int = 4,
    on_task_complete: Optional[Callable] = None,
) -> List[ComparisonResult]:
    """Run tasks via ProcessPoolExecutor (workers=1 for sequential fallback)."""
    if not tasks:
        return []

    workers = max(1, min(workers, len(tasks)))
    total = len(tasks)
    results: List[ComparisonResult] = []

    if workers == 1:
        for i, task in enumerate(tasks):
            result = process_single_task(task)
            results.append(result)
            if on_task_complete:
                on_task_complete(result, i + 1, total)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            future_map = {
                pool.submit(process_single_task, task): task
                for task in tasks
            }

            completed = 0
            for future in as_completed(future_map):
                completed += 1
                try:
                    result = future.result()
                except Exception as error:
                    task = future_map[future]
                    result = ComparisonResult(
                        task_id=task.task_id,
                        task_type=task.task_type,
                        success=False,
                        error=f"Process error: {error}",
                    )
                results.append(result)
                if on_task_complete:
                    on_task_complete(result, completed, total)

    return results
