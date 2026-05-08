"""Summary-plot orchestration: data prep, registry, output writer.

Loads the comparison CSV, aggregates per-peptidoform statistics, and delegates
to the per-plot-type modules (``summary_dual_axis``, ``summary_ridgeline``).
Shared constants (metric list, colour palette, plot-type registry) live in
``_summary_shared`` so both renderers stay in sync.
"""
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mspecviz.common.defaults import DEFAULTS
from mspecviz.plotting._summary_shared import METRIC_COLUMNS, _metric_label
from mspecviz.plotting.summary_dual_axis import plot_classic_dual_axis
from mspecviz.plotting.summary_ridgeline import plot_ridgeline

logger = logging.getLogger(__name__)

_PLOT_REGISTRY = {
    "dual_axis": plot_classic_dual_axis,
    "ridgeline": plot_ridgeline,
}


def prepare_summary_data(
    source: str | Path | pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and aggregate comparison data into per-peptidoform statistics.

    Returns (raw_frame, aggregated_frame).
    """
    if isinstance(source, (str, Path)):
        raw_frame = pd.read_csv(source, low_memory=False)
    else:
        raw_frame = source.copy()

    if metric not in raw_frame.columns:
        raise ValueError(
            f"Metric '{metric}' not found. Available: "
            f"{[c for c in raw_frame.columns if c in METRIC_COLUMNS]}"
        )

    # Build peptidoform label. ``Modification`` already stores the canonical
    # ``Peptidoform.modified_sequence`` (e.g. ``"PEPTM[Ox]IDEK"``) for modified
    # PSMs, or ``""``/``"(unmodified)"`` for unmodified — fall back to the bare
    # peptide in those cases.
    modification = raw_frame["Modification"].fillna("").astype(str)
    inline = modification.mask(
        modification.isin(["", "(unmodified)"]), raw_frame["Peptide"].astype(str),
    )
    raw_frame["peptidoform"] = inline + "/" + raw_frame["Charge"].astype(str)

    raw_frame = raw_frame.dropna(subset=[metric]).copy()

    grouped = raw_frame.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg([
        "count", "median", "mean", "std",
        ("q25", lambda x: x.quantile(0.25)),
        ("q75", lambda x: x.quantile(0.75)),
    ])
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["charge"] = grouped["Charge"].agg(lambda x: x.mode().iloc[0])

    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_frame, agg


def generate_all_summary_plots(
    source: str | Path | pd.DataFrame,
    output_dir: str | Path,
    metric: str = DEFAULTS.summary_metric,
    plot_types: list[str] | None = None,
    formats: list[str] | None = None,
    export_resolution: int = 300,
) -> dict[str, dict[str, str]]:
    """Generate all (or selected) summary plots and save them.

    Returns ``{plot_type: {format: file_path}}``.
    """
    if plot_types is None:
        plot_types = list(DEFAULTS.summary_plot_types)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    raw_frame, aggregated_frame = prepare_summary_data(source, metric=metric)
    logger.info(f"Summary data: {len(aggregated_frame)} peptidoforms, {len(raw_frame)} rows, metric={metric}")

    results: dict[str, dict[str, str]] = {}
    for plot_type in plot_types:
        func = _PLOT_REGISTRY.get(plot_type)
        if func is None:
            logger.warning("Unknown plot type: %s (skipping)", plot_type)
            continue

        try:
            fig = func(raw_frame, aggregated_frame, metric=metric,
                       export_resolution=export_resolution, title=None)
        except ValueError as error:
            logger.error("Failed to create %s plot: %s", plot_type, error)
            continue

        basename = ("ridgeline_r2_binned"
                    if plot_type == "ridgeline" and metric == DEFAULTS.summary_metric
                    else plot_type)
        paths: dict[str, str] = {}
        for output_format in formats:
            output_path = summary_dir / f"{basename}.{output_format}"
            try:
                fig.savefig(str(output_path), format=output_format, dpi=export_resolution,
                            bbox_inches="tight", facecolor="white")
                paths[output_format] = str(output_path)
            except OSError as error:
                logger.error("Failed to save %s.%s: %s", basename, output_format, error)
        plt.close(fig)
        if paths:
            results[plot_type] = paths
    return results
