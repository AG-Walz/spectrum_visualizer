"""Metric list and figure size helper shared by summary_dual_axis and summary_ridgeline 
used by core.plotting.__init__, summary_pipeline, summary_dual_axis and summary_ridgeline
"""
# Ordered comparison-metric columns expected in the summary DataFrame — consumed
# by summary_pipeline to select metric columns and re-exported via
# core.plotting.__init__. Order defines the layout order of summary plots.
METRIC_COLUMNS = [
    "R² match", "R² binned", "Cosine", "Entropy",
    "Dice", "Jaccard", "SharedRatio",
]

# Maps each raw METRIC_COLUMNS name to a verbose axis-label string used when
# rendering summary plots.
_METRIC_LABELS = {
    "R² match": "R² (match)",
    "R² binned": "R² (binned)",
    "Cosine": "Cosine similarity",
    "Entropy": "Spectral entropy similarity",
    "Dice": "Dice coefficient",
    "Jaccard": "Jaccard index",
    "SharedRatio": "Shared peak ratio",
}


def _metric_label(metric: str) -> str:
    """Return the human-readable axis label for metric.
    Falls back to the raw metric name itself when it is not a known key in
    _METRIC_LABELS, so an unrecognised metric still renders sensibly.
    """
    return _METRIC_LABELS.get(metric, metric)


def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    """Compute a figure height (in inches) that scales with the number of rows.
    Allocates per_row inches per variant row but never returns less than
    minimum, keeping small plots legible while letting tall ones grow.
    """
    return max(minimum, per_row * n)
