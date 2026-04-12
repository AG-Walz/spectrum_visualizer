"""
Parallel plot generation via ProcessPoolExecutor.

Main process annotates spectra, worker pool runs process_single_task()
(metrics + matplotlib/plotly), main process collects ComparisonResults.
"""
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.mirror_plot import (
    create_mirror_plot, create_single_spectrum_plot,
    create_plotly_mirror_plot, create_plotly_single_plot,
)
from src.spectrum_metrics import compute_all_metrics, format_metrics_text

logger = logging.getLogger(__name__)


@dataclass
class ComparisonTask:
    """All data needed to produce one plot in a worker process."""

    task_id: int
    task_type: str  # "mirror" or "single"

    # Annotated spectrum DataFrames (picklable)
    query_df: pd.DataFrame
    ref_df: Optional[pd.DataFrame] = None

    # Labels for the plot
    query_label: str = ""
    ref_label: str = ""
    plot_title: str = ""

    # Export paths (strings -- workers create subdirs as needed)
    folder_name: str = ""
    safe_filename: str = ""
    png_dir: str = ""
    svg_dir: str = ""
    pdf_dir: str = ""

    # Plot/export parameters
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    figsize: Tuple[int, int] = (12, 8)
    dpi: int = 300
    normalize: bool = True

    # Output path for interactive Plotly HTML (empty string = skip HTML output)
    html_plot_path: str = ""

    # Skip plot generation (metrics only -- plot already exists on disk)
    skip_plot: bool = False

    # Summary metadata (for CSV row)
    query_file: str = ""
    ref_file: str = ""
    peptide: str = ""
    modification: str = ""
    charge: int = 0

    def outputs_exist(self) -> bool:
        """Return True if every requested output format already exists on disk."""
        fmt_dirs = {"png": self.png_dir, "svg": self.svg_dir, "pdf": self.pdf_dir}
        for fmt, fmt_dir in fmt_dirs.items():
            if fmt in self.export_formats and fmt_dir:
                out = Path(fmt_dir) / self.folder_name / f"{self.safe_filename}.{fmt}"
                if not out.exists():
                    return False
        if "html" in self.export_formats and self.html_plot_path:
            if not Path(self.html_plot_path).exists():
                return False
        return True


@dataclass
class ComparisonResult:
    """Output from one worker invocation."""

    task_id: int
    task_type: str
    success: bool
    error: str = ""

    png_path: str = ""
    svg_path: str = ""
    pdf_path: str = ""
    html_path: str = ""

    # Dict with metrics for summary CSV (only for mirror tasks)
    summary_row: Optional[Dict] = None


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
                metrics_text = format_metrics_text(metrics)
                fig = create_mirror_plot(
                    upper_spectrum=task.query_df,
                    lower_spectrum=task.ref_df,
                    upper_label=task.query_label,
                    lower_label=task.ref_label,
                    title=task.plot_title,
                    metrics_text=metrics_text,
                    figsize=task.figsize,
                    dpi=task.dpi,
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
                figsize=task.figsize,
                dpi=task.dpi,
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
            save_kwargs["dpi"] = task.dpi
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
