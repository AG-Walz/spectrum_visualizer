"""
Parallel plot generation for MSpecViz.

Moves the heavy per-plot work (metrics computation, matplotlib figure creation,
savefig) into worker processes via ProcessPoolExecutor.  The annotation step
stays in the main process (fast, ~10 ms) while the expensive steps (~400 ms
each) run in parallel.

Architecture:
    main process:  stream queries → annotate → build ComparisonTask list
    worker pool:   process_single_task(task) → ComparisonResult
    main process:  collect results, update progress, write summary
"""
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


# ============================================================================
# Data containers — must be picklable for multiprocessing
# ============================================================================

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

    # Export paths (strings — workers create subdirs as needed)
    folder_name: str = ""
    safe_filename: str = ""
    png_dir: str = ""
    svg_dir: str = ""

    # Plot/export parameters
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    figsize: Tuple[int, int] = (12, 8)
    dpi: int = 300
    normalize: bool = True

    # Skip plot generation (metrics only — plot already exists on disk)
    skip_plot: bool = False

    # Summary metadata (for CSV row)
    query_file: str = ""
    ref_file: str = ""
    peptide: str = ""
    modification: str = ""
    charge: int = 0


@dataclass
class ComparisonResult:
    """Output from one worker invocation."""

    task_id: int
    task_type: str
    success: bool
    error: str = ""

    png_path: str = ""
    svg_path: str = ""

    # Dict with metrics for summary CSV (only for mirror tasks)
    summary_row: Optional[Dict] = None


# ============================================================================
# Worker function — runs in a child process
# ============================================================================

def process_single_task(task: ComparisonTask) -> ComparisonResult:
    """
    Top-level worker: compute metrics, create figure, save to disk.

    If task.skip_plot is True, only computes metrics (for summary CSV)
    without creating or saving a figure — used when the plot file already
    exists on disk from a previous run.

    Imports are done inside the function so each worker process initialises
    its own matplotlib backend and module references.
    """
    try:
        summary_row = None
        png_path = ""
        svg_path = ""

        if task.task_type == "mirror" and task.ref_df is not None:
            # --- Metrics (always computed for summary CSV) ---
            from src.spectrum_metrics import compute_all_metrics

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
                # --- Full mirror plot ---
                import matplotlib
                matplotlib.use("Agg")
                from src.mirror_plot import create_mirror_plot
                from src.spectrum_metrics import format_metrics_text

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
                png_path, svg_path = _save_task_figure(fig, task)
                import matplotlib.pyplot as plt
                plt.close(fig)
                del fig

        elif not task.skip_plot:
            # --- Single spectrum plot (no metrics) ---
            import matplotlib
            matplotlib.use("Agg")
            from src.mirror_plot import create_single_spectrum_plot

            fig = create_single_spectrum_plot(
                spectrum=task.query_df,
                title=task.plot_title,
                label=task.query_label,
                figsize=task.figsize,
                dpi=task.dpi,
                normalize=task.normalize,
            )
            png_path, svg_path = _save_task_figure(fig, task)
            import matplotlib.pyplot as plt
            plt.close(fig)
            del fig

        return ComparisonResult(
            task_id=task.task_id,
            task_type=task.task_type,
            success=True,
            png_path=png_path,
            svg_path=svg_path,
            summary_row=summary_row,
        )

    except Exception as e:
        return ComparisonResult(
            task_id=task.task_id,
            task_type=task.task_type,
            success=False,
            error=str(e),
        )


def _save_task_figure(fig, task: ComparisonTask) -> Tuple[str, str]:
    """Save a matplotlib figure to PNG and/or SVG. Returns (png_path, svg_path)."""
    png_path = ""
    svg_path = ""

    if "png" in task.export_formats and task.png_dir:
        variant_dir = Path(task.png_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        out = variant_dir / f"{task.safe_filename}.png"
        fig.savefig(str(out), dpi=task.dpi, bbox_inches="tight", facecolor="white")
        png_path = str(out)

    if "svg" in task.export_formats and task.svg_dir:
        variant_dir = Path(task.svg_dir) / task.folder_name
        variant_dir.mkdir(parents=True, exist_ok=True)
        out = variant_dir / f"{task.safe_filename}.svg"
        fig.savefig(str(out), format="svg", bbox_inches="tight", facecolor="white")
        svg_path = str(out)

    return png_path, svg_path


# ============================================================================
# Orchestrator — called from the main process
# ============================================================================

def run_parallel_comparisons(
    tasks: List[ComparisonTask],
    workers: int = 4,
    on_task_complete: Optional[Callable] = None,
) -> List[ComparisonResult]:
    """
    Process a list of ComparisonTasks, optionally in parallel.

    Args:
        tasks: List of tasks to process.
        workers: Number of worker processes.  Clamped to len(tasks).
                 If 1, runs sequentially (useful for debugging).
        on_task_complete: Optional callback(result, completed, total)
                         called from the main thread after each task finishes.

    Returns:
        List of ComparisonResult in completion order.
    """
    if not tasks:
        return []

    workers = max(1, min(workers, len(tasks)))
    total = len(tasks)
    results: List[ComparisonResult] = []

    if workers == 1:
        # Sequential fallback — no subprocess overhead, easier to debug
        for i, task in enumerate(tasks):
            result = process_single_task(task)
            results.append(result)
            if on_task_complete:
                on_task_complete(result, i + 1, total)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

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
