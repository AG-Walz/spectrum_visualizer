"""Summary-plot orchestration: data prep, registry, output writer.

Loads the comparison CSV, aggregates per-peptidoform statistics, and delegates
to the per-plot-type modules (``summary_dual_axis``, ``summary_ridgeline``).
Shared constants (metric list, colour palette, plot-type registry) live in
``_summary_shared`` so both renderers stay in sync.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mspecviz.plotting._summary_shared import (
    DEFAULT_METRIC,
    MAX_SCATTER_PER_GROUP,
    METRIC_COLUMNS,
    PLOT_TYPES,
    SUMMARY_COLORS,
    _auto_height,
    _metric_label,
    _style_ax,
    _subsample_group,
)
from mspecviz.plotting.summary_dual_axis import plot_classic_dual_axis
from mspecviz.plotting.summary_ridgeline import plot_ridgeline
from mspecviz.search.models import format_inline_sequence

logger = logging.getLogger(__name__)

_PLOT_REGISTRY = {
    "dual_axis": plot_classic_dual_axis,
    "ridgeline": plot_ridgeline,
}


def prepare_summary_data(
    source: Union[str, Path, pd.DataFrame],
    metric: str = DEFAULT_METRIC,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
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

    # Build peptidoform label
    modification = raw_frame["Modification"].fillna("").astype(str)
    raw_frame["peptidoform"] = (
        raw_frame.apply(
            lambda row: format_inline_sequence(row["Peptide"], modification[row.name]),
            axis=1,
        )
        + "/"
        + raw_frame["Charge"].astype(str)
    )

    raw_frame = raw_frame.dropna(subset=[metric]).copy()

    grouped = raw_frame.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg(["count", "median", "mean", "std"])
    agg.columns = ["count", "median", "mean", "std"]
    agg["q25"] = grouped[metric].quantile(0.25)
    agg["q75"] = grouped[metric].quantile(0.75)
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])

    charge_mode = grouped["Charge"].agg(lambda x: x.mode().iloc[0] if len(x) > 0 else 0)
    agg["charge"] = charge_mode

    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_frame, agg


def generate_all_summary_plots(
    source: Union[str, Path, pd.DataFrame],
    output_dir: Union[str, Path],
    metric: str = DEFAULT_METRIC,
    plot_types: Optional[List[str]] = None,
    formats: Optional[List[str]] = None,
    export_resolution: int = 300,
    title_prefix: str = "",
) -> Dict[str, Dict[str, str]]:
    """Generate all (or selected) summary plots and save them.

    Returns ``{plot_type: {format: file_path}}``.
    """
    if plot_types is None:
        plot_types = list(PLOT_TYPES)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    raw_frame, aggregated_frame = prepare_summary_data(source, metric=metric)
    logger.info("Summary data: %d peptidoforms, %d rows, metric=%s",
                len(aggregated_frame), len(raw_frame), metric)

    results: Dict[str, Dict[str, str]] = {}

    for pt in plot_types:
        func = _PLOT_REGISTRY.get(pt)
        if func is None:
            logger.warning("Unknown plot type: %s (skipping)", pt)
            continue

        plot_title = (f"{title_prefix} -- {_metric_label(metric)}".strip(" -")
                      if title_prefix else None)

        try:
            fig = func(raw_frame, aggregated_frame, metric=metric, export_resolution=export_resolution, title=plot_title)
        except ValueError as exc:
            logger.error("Failed to create %s plot: %s", pt, exc)
            continue

        paths: Dict[str, str] = {}
        filename = f"ridgeline_r2_binned" if pt == "ridgeline" and metric == DEFAULT_METRIC else pt
        for fmt in formats:
            out_path = summary_dir / f"{filename}.{fmt}"
            try:
                fig.savefig(str(out_path), format=fmt, dpi=export_resolution,
                            bbox_inches="tight", facecolor="white")
                paths[fmt] = str(out_path)
            except OSError as exc:
                logger.error("Failed to save %s.%s: %s", pt, fmt, exc)

        plt.close(fig)
        if paths:
            results[pt] = paths

    return results
