"""Constants + helpers shared by the summary renderers.

Lives in its own module so ``summary_dual_axis`` and ``summary_ridgeline`` can
import without creating a circular dependency back to ``summary_pipeline``
(which imports the renderers at module top).
"""
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
