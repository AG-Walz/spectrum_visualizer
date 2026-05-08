"""Set-based + entropy similarity metrics, full-metric aggregator, and formatting.

Set metrics (Dice / Jaccard / shared-ratio) operate on the *set* of shared
annotations rather than their intensities. This file owns the
``compute_all_metrics`` aggregator that the rest of the pipeline calls — it
combines intensity metrics from :mod:`metrics_intensity` with the set
metrics here, and produces the dict consumed by ``format_metrics_text``.
"""
import logging
import math

import numpy as np
import pandas as pd

from mspecviz.annotation.metrics_intensity import (
    _cosine_from_arrays,
    _get_matching_intensities,
    _r2_from_arrays,
    compute_r2_binned,
)

try:
    import ms_entropy
except ImportError:
    ms_entropy = None

logger = logging.getLogger(__name__)


def compute_spectral_entropy_similarity(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> float:
    """Spectral entropy similarity (Li et al., Nature Methods 2021)."""
    if query_frame.empty or reference_frame.empty:
        return float("nan")

    if ms_entropy is None:
        logger.warning("ms_entropy not installed -- entropy similarity unavailable")
        return float("nan")

    q_spec = np.column_stack([query_frame["mz"].values, query_frame["intensity"].values])
    r_spec = np.column_stack([reference_frame["mz"].values, reference_frame["intensity"].values])

    try:
        similarity = ms_entropy.calculate_entropy_similarity(
            q_spec, r_spec, ms2_tolerance_in_da=0.02
        )
        return float(similarity)
    except (ValueError, FloatingPointError) as error:
        logger.warning(f"Entropy similarity failed: {error}")
        return float("nan")


def compute_all_metrics(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> dict[str, float | int]:
    """Compute all 7 comparison metrics plus context counts.

    Returns dict with keys:
        Metrics:  R2m, R2b, Cos, Ent, Dice, Jacc, ShR
        Counts:   total_q, total_r, ann_q, ann_r, shared
    """
    query_intensities, reference_intensities = _get_matching_intensities(query_frame, reference_frame)
    n_shared = len(query_intensities)

    ann_q = int(query_frame["is_annotated"].sum())
    ann_r = int(reference_frame["is_annotated"].sum())

    denom_dice = ann_q + ann_r
    denom_jacc = ann_q + ann_r - n_shared
    denom_shr = min(ann_q, ann_r)

    return {
        "R2m": _r2_from_arrays(query_intensities, reference_intensities),
        "R2b": compute_r2_binned(query_frame, reference_frame),
        "Cos": _cosine_from_arrays(query_intensities, reference_intensities),
        "Ent": compute_spectral_entropy_similarity(query_frame, reference_frame),
        "Dice": float(2 * n_shared / denom_dice) if denom_dice else float("nan"),
        "Jacc": float(n_shared / denom_jacc) if denom_jacc else float("nan"),
        "ShR": float(n_shared / denom_shr) if denom_shr else float("nan"),
        "total_q": len(query_frame),
        "total_r": len(reference_frame),
        "ann_q": ann_q,
        "ann_r": ann_r,
        "shared": n_shared,
    }


def format_metrics_text(metrics: dict[str, float | int]) -> str:
    """Format metrics as multi-line string for plot overlay."""
    def _fmt(value) -> str:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "N/A"
        return f"{value:.3f}"

    return (
        f"R² Match: {_fmt(metrics.get('R2m'))}\n"
        f"R² Binned: {_fmt(metrics.get('R2b'))}\n"
        f"Cosine: {_fmt(metrics.get('Cos'))}\n"
        f"Entropy: {_fmt(metrics.get('Ent'))}\n"
        f"Dice: {_fmt(metrics.get('Dice'))}\n"
        f"Jaccard: {_fmt(metrics.get('Jacc'))}\n"
        f"Shared: {_fmt(metrics.get('ShR'))}"
    )


# Mapping from short keys to display column names (used by summary export)
METRIC_SHORT_TO_DISPLAY = {
    "R2m": "R² match",
    "R2b": "R² binned",
    "Cos": "Cosine",
    "Ent": "Entropy",
    "Dice": "Dice",
    "Jacc": "Jaccard",
    "ShR": "SharedRatio",
}
