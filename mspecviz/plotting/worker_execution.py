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


def process_single_task(task: ComparisonTask) -> ComparisonResult:
    """Worker: compute metrics, create figure, save to disk. Skips plot if skip_plot is set."""
    try:
        summary_row = None
        png_path = ""
        svg_path = ""
        pdf_path = ""
        html_path = ""

        if task.task_type == "mirror" and task.ref_df is not None:
            # --- Metrics (always computed for summary CSV) ---
            metrics = compute_all_metrics(task.query_df, task.ref_df)

            summary_row = {
                "query_file": task.query_file,
                "ref_file": task.ref_file,
                "peptide": task.peptide,
                "modification": task.modification,
                "charge": task.charge,
                **metrics,
            }

            if not task.skip_plot:
                # --- Full mirror plot (Matplotlib, batch export) ---
                metrics_text = task.metrics_text or format_metrics_text(metrics)
                fig = create_mirror_plot(
                    upper_spectrum=task.query_df,
                    lower_spectrum=task.ref_df,
                    upper_label=task.query_label,
                    lower_label=task.ref_label,
                    title=task.plot_title,
                    metrics_text=metrics_text,
                    figure_size=task.figure_size,
                    export_dpi=task.export_dpi,
                    normalize=task.normalize,
                )
                png_path, svg_path, pdf_path = _save_task_figure(fig, task)
                plt.close(fig)
                del fig

            # --- Interactive HTML (Plotly, optional) ---
            if task.html_plot_path:
                html_path = _save_task_html(task)

        elif not task.skip_plot:
            # --- Single spectrum plot (no metrics) ---
            fig = create_single_spectrum_plot(
                spectrum=task.query_df,
                title=task.plot_title,
                label=task.query_label,
                figure_size=task.figure_size,
                export_dpi=task.export_dpi,
                normalize=task.normalize,
            )
            png_path, svg_path, pdf_path = _save_task_figure(fig, task)
            plt.close(fig)
            del fig

            # --- Interactive HTML for single-spectrum plots ---
            if task.html_plot_path:
                html_path = _save_task_html(task)

        return ComparisonResult(
            task_id=task.task_id,
            task_type=task.task_type,
            success=True,
            png_path=png_path,
            svg_path=svg_path,
            pdf_path=pdf_path,
            html_path=html_path,
            summary_row=summary_row,
        )

    except Exception as e:
        return ComparisonResult(
            task_id=task.task_id,
            task_type=task.task_type,
            success=False,
            error=str(e),
        )


# Format -> savefig keyword and extension mapping
_FORMAT_SAVE_KWARGS: Dict[str, dict] = {
    "png": {"format": "png"},
    "svg": {"format": "svg"},
    "pdf": {"format": "pdf"},
}


def _save_task_figure(fig, task: ComparisonTask) -> Tuple[str, str, str]:
    """Save a matplotlib figure to each requested format. Returns (png_path, svg_path, pdf_path)."""
    fmt_dirs = {"png": task.png_dir, "svg": task.svg_dir, "pdf": task.pdf_dir}
    paths = {"png": "", "svg": "", "pdf": ""}

    for fmt, fmt_dir in fmt_dirs.items():
        if fmt not in task.export_formats or not fmt_dir:
            continue
        variant_dir = Path(fmt_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        out = variant_dir / f"{task.safe_filename}.{fmt}"
        save_kwargs = {"bbox_inches": "tight", "facecolor": "white", **_FORMAT_SAVE_KWARGS[fmt]}
        if fmt == "png":
            save_kwargs["dpi"] = task.export_dpi
        fig.savefig(str(out), **save_kwargs)
        paths[fmt] = str(out)

    return paths["png"], paths["svg"], paths["pdf"]


def _save_task_html(task: ComparisonTask) -> str:
    """Save interactive Plotly HTML plot. Returns path or "" on failure."""
    if not task.html_plot_path:
        return ""

    try:
        out = Path(task.html_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)

        if task.task_type == "mirror" and task.ref_df is not None:
            fig = create_plotly_mirror_plot(
                upper_spectrum=task.query_df,
                lower_spectrum=task.ref_df,
                title=task.plot_title,
                normalize=task.normalize,
            )
        else:
            fig = create_plotly_single_plot(
                spectrum=task.query_df,
                title=task.plot_title,
                normalize=task.normalize,
            )

        fig.write_html(str(out), include_plotlyjs="cdn")
        return str(out)

    except Exception as e:
        logger.warning(f"HTML plot save failed for task {task.task_id}: {e}")
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
                except Exception as e:
                    task = future_map[future]
                    result = ComparisonResult(
                        task_id=task.task_id,
                        task_type=task.task_type,
                        success=False,
                        error=f"Process error: {e}",
                    )
                results.append(result)
                if on_task_complete:
                    on_task_complete(result, completed, total)

    return results
