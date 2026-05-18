"""Metric list, pretty labels, and height helper shared by ``summary_dual_axis`` and ``summary_ridgeline`` (kept here to avoid a circular import back through ``summary_pipeline``)."""
METRIC_COLUMNS = [
    "R² match", "R² binned", "Cosine", "Entropy",
    "Dice", "Jaccard", "SharedRatio",
]

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
    return _METRIC_LABELS.get(metric, metric)


def _auto_height(n: int, per_row: float = 0.45, minimum: float = 5.0) -> float:
    return max(minimum, per_row * n)
