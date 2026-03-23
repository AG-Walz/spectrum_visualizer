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
    pdf_dir: str = ""

    # Plot/export parameters
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    figsize: Tuple[int, int] = (12, 8)
    dpi: int = 300
    normalize: bool = True

    # Output path for interactive Plotly HTML (empty string = skip HTML output)
    html_plot_path: str = ""

    # Skip plot generation (metrics only — plot already exists on disk)
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


# Keys written by SpectrumAnnotator.addIonMatchStatistics() as PeptideHit MetaValues
_IMS_STAT_KEYS = [
    "NTermIonCurrentRatio", "CTermIonCurrentRatio", "TIC", "MIC",
    "annotated_tic_fraction", "matched_prefix_ions_fraction",
    "matched_suffix_ions_fraction", "max_series_type", "max_series_charge",
    "max_series_size", "coverage", "annotated_ion_fraction", "ms2_scan_index",
    "precursor_mz_error_ppm", "fragment_mass_error_mean_ppm",
    "fragment_mass_error_stdev_ppm", "median_fragment_error_ppm",
]


def compute_ion_match_statistics(
    annotated_df: pd.DataFrame,
    peptidoform_str: str,
    precursor_charge: int,
    precursor_mz: float = 0.0,
) -> dict:
    """
    Compute 17 SpectrumAnnotator ion-match statistics as supplementary metrics.

    Wraps oms.SpectrumAnnotator.addIonMatchStatistics().  Results are returned
    with an ``ims_`` prefix to avoid collision with existing metric columns.

    This is called AFTER the primary annotation pipeline (annotate_spectrum /
    prepare_spectrum_for_plot) and does NOT replace it.  The primary pipeline
    uses SpectrumAlignment + manual pairs loop which correctly handles charge ≥3;
    this function only produces supplementary statistics.

    Args:
        annotated_df: Annotated spectrum DataFrame (mz, intensity columns required)
        peptidoform_str: ProForma peptidoform string
        precursor_charge: Precursor charge state
        precursor_mz: Precursor m/z (used for precursor_mz_error_ppm stat)

    Returns:
        Dict of {ims_<key>: value} for all 17 stats; None values on failure.
    """
    try:
        import pyopenms as oms
        from src.peak_annotation import peptidoform_to_aasequence

        # Build MSSpectrum from annotated DataFrame
        spec = oms.MSSpectrum()
        spec.set_peaks([
            annotated_df["mz"].tolist(),
            annotated_df["intensity"].tolist(),
        ])
        spec.sortByPosition()
        if precursor_mz > 0:
            precursor = oms.Precursor()
            precursor.setMZ(precursor_mz)
            precursor.setCharge(precursor_charge)
            spec.setPrecursors([precursor])

        # Build PeptideIdentification with PeptideHit
        aa_seq = peptidoform_to_aasequence(peptidoform_str)
        ph = oms.PeptideHit()
        ph.setSequence(aa_seq)
        ph.setCharge(precursor_charge)
        pi = oms.PeptideIdentification()
        pi.setHits([ph])

        # Configure a minimal TSG and SpectrumAlignment for addIonMatchStatistics
        tsg = oms.TheoreticalSpectrumGenerator()
        tsg_params = tsg.getParameters()
        tsg_params.setValue("add_b_ions", "true")
        tsg_params.setValue("add_y_ions", "true")
        tsg_params.setValue("add_metainfo", "true")
        tsg.setParameters(tsg_params)

        sa = oms.SpectrumAlignment()
        sa_params = sa.getParameters()
        sa_params.setValue("tolerance", 0.02)
        sa_params.setValue("is_relative_tolerance", "false")
        sa.setParameters(sa_params)

        annotator = oms.SpectrumAnnotator()
        # Signature: addIonMatchStatistics(pi, spec, tg, sa)
        annotator.addIonMatchStatistics(pi, spec, tsg, sa)

        updated_ph = pi.getHits()[0]
        stats = {}
        for key in _IMS_STAT_KEYS:
            try:
                stats[f"ims_{key}"] = updated_ph.getMetaValue(key)
            except Exception:
                stats[f"ims_{key}"] = None
        return stats

    except Exception as e:
        logger.warning(f"compute_ion_match_statistics failed for '{peptidoform_str}': {e}")
        return {f"ims_{k}": None for k in _IMS_STAT_KEYS}


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
        pdf_path = ""
        html_path = ""

        if task.task_type == "mirror" and task.ref_df is not None:
            # --- Metrics (always computed for summary CSV) ---
            from src.spectrum_metrics import compute_all_metrics

            metrics = compute_all_metrics(task.query_df, task.ref_df)

            # Supplementary ims_ stats via SpectrumAnnotator.addIonMatchStatistics()
            ims_stats = compute_ion_match_statistics(
                task.query_df, task.peptide, task.charge,
            )

            summary_row = {
                "query_file": task.query_file,
                "ref_file": task.ref_file,
                "peptide": task.peptide,
                "modification": task.modification,
                "charge": task.charge,
                **metrics,
                **ims_stats,
            }

            if not task.skip_plot:
                # --- Full mirror plot (Matplotlib, batch export) ---
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
                png_path, svg_path, pdf_path = _save_task_figure(fig, task)
                import matplotlib.pyplot as plt
                plt.close(fig)
                del fig

            # --- Interactive HTML (pyopenms-viz Plotly, optional) ---
            if task.html_plot_path:
                html_path = _save_task_html(task)

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
            png_path, svg_path, pdf_path = _save_task_figure(fig, task)
            import matplotlib.pyplot as plt
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


# Format → savefig keyword and extension mapping
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
    """
    Save an interactive Plotly HTML plot alongside PNG/SVG outputs.

    Uses pyopenms-viz ms_plotly backend.  The HTML uses include_plotlyjs='cdn'
    (~10 KB file size) so an internet connection is required to view it.

    Returns the saved path as a string, or "" on failure.
    """
    if not task.html_plot_path:
        return ""

    try:
        out = Path(task.html_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)

        if task.task_type == "mirror" and task.ref_df is not None:
            from src.mirror_plot import create_interactive_mirror_plot
            fig = create_interactive_mirror_plot(
                upper_spectrum=task.query_df,
                lower_spectrum=task.ref_df,
                title=task.plot_title,
                relative_intensity=task.normalize,
            )
        else:
            from src.mirror_plot import create_interactive_single_plot
            fig = create_interactive_single_plot(
                spectrum=task.query_df,
                title=task.plot_title,
                relative_intensity=task.normalize,
            )

        fig.write_html(str(out), include_plotlyjs="cdn")
        return str(out)

    except Exception as e:
        logger.warning(f"HTML plot save failed for task {task.task_id}: {e}")
        return ""


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
