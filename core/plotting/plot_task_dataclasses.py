"""``ComparisonTask`` / ``ComparisonResult`` dataclasses passed between the orchestrator and the parallel render workers (constructed by ``core.export.export_pipeline`` / ``cli.main``, consumed by ``plot_worker``)."""
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ComparisonTask:
    """Picklable bundle of everything a worker needs to render one plot (spectra, labels, output paths, formats)."""

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

    # Pre-formatted metrics box: lets the caller compute metrics once in the
    # main process and reuse the text across workers instead of recomputing.
    metrics_text: str = ""

    # Summary metadata (for CSV row)
    query_file: str = ""
    reference_file: str = ""
    peptide: str = ""
    modification: str = ""
    charge: int = 0


@dataclass
class ComparisonResult:
    """Worker output: success flag, on-disk paths for each rendered format, and the summary-CSV row for mirror tasks."""

    task_id: int
    task_type: str
    success: bool
    error: str = ""

    png_path: str = ""
    svg_path: str = ""
    pdf_path: str = ""
    html_path: str = ""

    # Per-task metrics row for the summary CSV (mirror tasks only)
    summary_row: Optional[dict] = None
