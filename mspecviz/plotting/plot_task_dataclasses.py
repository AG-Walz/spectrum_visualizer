"""Task / result dataclasses for parallel plot generation.

Kept separate from ``worker_execution`` so they remain lightweight to import
(no matplotlib / plotly dependency) for callers that only need to construct
tasks — e.g. the CLI builder in ``mspecviz.pipeline``.
"""
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ComparisonTask:
    """All data needed to produce one plot in a worker process."""

    task_id: int
    task_type: str  # "mirror" or "single"

    # Annotated spectrum DataFrames (picklable)
    query_frame: pd.DataFrame
    reference_frame: Optional[pd.DataFrame] = None

    # Labels for the plot
    query_label: str = ""
    reference_label: str = ""
    plot_title: str = ""

    # Export paths (strings -- workers create subdirs as needed)
    folder_name: str = ""
    safe_filename: str = ""
    png_dir: str = ""
    svg_dir: str = ""
    pdf_dir: str = ""

    # Plot/export parameters
    export_formats: list[str] = field(default_factory=lambda: ["png"])
    figure_size: tuple[int, int] = (12, 8)
    export_resolution: int = 300

    # Output path for interactive Plotly HTML (empty string = skip HTML output)
    html_plot_path: str = ""

    # Skip plot generation (metrics only -- plot already exists on disk)
    skip_plot: bool = False

    # Pre-formatted metrics box text. When non-empty, the worker uses this
    # directly instead of recomputing via ``format_metrics_text`` — lets the
    # caller (e.g. the CLI) compute metrics once in the main process and
    # avoid re-doing the work per worker.
    metrics_text: str = ""

    # Summary metadata (for CSV row)
    query_file: str = ""
    reference_file: str = ""
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
    pdf_path: str = ""
    html_path: str = ""

    # Dict with metrics for summary CSV (only for mirror tasks)
    summary_row: Optional[dict] = None
