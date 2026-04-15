"""Summary-plot orchestration: data prep, registry, output writer.

Loads the comparison CSV, aggregates per-peptidoform statistics, and delegates
to the per-plot-type modules (``summary_dual_axis``, ``summary_ridgeline``).
Shared constants (metric list, colour palette, plot-type registry) live here
so both renderers stay in sync.
"""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mspecviz.search.models import format_inline_sequence

logger = logging.getLogger(__name__)


METRIC_COLUMNS = [
    "R\u00b2 match",
    "R\u00b2 binned",
    "Cosine",
    "Entropy",
    "Dice",
    "Jaccard",
    "SharedRatio",
]

DEFAULT_METRIC = "R\u00b2 binned"

SUMMARY_COLORS = {
    "bar": "#E8945A",
    "bar_edge": "#D07030",
    "box": "#2C3E50",
    "box_fill": "#BDC3C7",
    "scatter": "#2C3E50",
}

PLOT_TYPES = ["dual_axis", "ridgeline"]

MAX_SCATTER_PER_GROUP = 200


def prepare_summary_data(
    source: Union[str, Path, pd.DataFrame],
    metric: str = DEFAULT_METRIC,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and aggregate comparison data into per-peptidoform statistics.

    Returns (raw_df, agg_df).
    """
    if isinstance(source, (str, Path)):
        raw_df = pd.read_csv(source, low_memory=False)
    else:
        raw_df = source.copy()

    if metric not in raw_df.columns:
        raise ValueError(
            f"Metric '{metric}' not found. Available: "
            f"{[c for c in raw_df.columns if c in METRIC_COLUMNS]}"
        )

    # Build peptidoform label
    mod = raw_df["Modification"].fillna("").astype(str)
    raw_df["peptidoform"] = (
        raw_df.apply(
            lambda row: format_inline_sequence(row["Peptide"], mod[row.name]),
            axis=1,
        )
        + "/"
        + raw_df["Charge"].astype(str)
    )

    raw_df = raw_df.dropna(subset=[metric]).copy()

    grouped = raw_df.groupby("peptidoform", sort=False)
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
    return raw_df, agg


def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    return max(minimum, per_row * n)


def _style_ax(ax, remove_top=True, remove_right=True):
    if remove_top:
        ax.spines["top"].set_visible(False)
    if remove_right:
        ax.spines["right"].set_visible(False)


def _subsample_group(group_df, metric, max_n=MAX_SCATTER_PER_GROUP, seed=42):
    vals = group_df[metric].values
    if len(vals) <= max_n:
        return vals
    rng = np.random.RandomState(seed)
    return rng.choice(vals, size=max_n, replace=False)


def _metric_label(metric: str) -> str:
    labels = {
        "R\u00b2 match": "R\u00b2 (match)",
        "R\u00b2 binned": "R\u00b2 (binned)",
        "Cosine": "Cosine similarity",
        "Entropy": "Spectral entropy similarity",
        "Dice": "Dice coefficient",
        "Jaccard": "Jaccard index",
        "SharedRatio": "Shared peak ratio",
    }
    return labels.get(metric, metric)


def generate_all_summary_plots(
    source: Union[str, Path, pd.DataFrame],
    output_dir: Union[str, Path],
    metric: str = DEFAULT_METRIC,
    plot_types: Optional[List[str]] = None,
    formats: Optional[List[str]] = None,
    dpi: int = 300,
    title_prefix: str = "",
) -> Dict[str, Dict[str, str]]:
    """Generate all (or selected) summary plots and save them.

    Returns ``{plot_type: {format: file_path}}``.
    """
    # Imported here to avoid a circular import — the per-plot modules depend
    # on helpers defined above.
    from mspecviz.plotting.summary_dual_axis import plot_classic_dual_axis
    from mspecviz.plotting.summary_ridgeline import plot_ridgeline

    _PLOT_REGISTRY = {
        "dual_axis": plot_classic_dual_axis,
        "ridgeline": plot_ridgeline,
    }

    if plot_types is None:
        plot_types = list(PLOT_TYPES)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    raw_df, agg_df = prepare_summary_data(source, metric=metric)
    logger.info("Summary data: %d peptidoforms, %d rows, metric=%s",
                len(agg_df), len(raw_df), metric)

    results: Dict[str, Dict[str, str]] = {}

    for pt in plot_types:
        func = _PLOT_REGISTRY.get(pt)
        if func is None:
            logger.warning("Unknown plot type: %s (skipping)", pt)
            continue

        plot_title = (f"{title_prefix} -- {_metric_label(metric)}".strip(" -")
                      if title_prefix else None)

        try:
            fig = func(raw_df, agg_df, metric=metric, dpi=dpi, title=plot_title)
        except Exception as exc:
            logger.error("Failed to create %s plot: %s", pt, exc)
            continue

        paths: Dict[str, str] = {}
        filename = f"ridgeline_r2_binned" if pt == "ridgeline" and metric == DEFAULT_METRIC else pt
        for fmt in formats:
            out_path = summary_dir / f"{filename}.{fmt}"
            try:
                fig.savefig(str(out_path), format=fmt, dpi=dpi,
                            bbox_inches="tight", facecolor="white")
                paths[fmt] = str(out_path)
            except Exception as exc:
                logger.error("Failed to save %s.%s: %s", pt, fmt, exc)

        plt.close(fig)
        if paths:
            results[pt] = paths

    return results
