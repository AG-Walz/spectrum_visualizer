"""Shared constants and helper functions for summary plot renderers.

Extracted here so that ``summary_dual_axis`` and ``summary_ridgeline`` can
import from this module without creating a circular dependency with
``summary_orchestration``.
"""
import numpy as np

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
