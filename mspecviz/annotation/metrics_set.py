"""Set-based + entropy similarity metrics, full-metric aggregator, and formatting.

Set metrics (Dice / Jaccard / shared-ratio) operate on the *set* of shared
annotations rather than their intensities. This file also owns the
``compute_all_metrics`` aggregator that the rest of the pipeline calls — it
combines intensity metrics from :mod:`metrics_intensity` with the set
metrics here, and produces the dict consumed by ``format_metrics_text``.
"""
import logging
import math
from pathlib import Path
from typing import Dict, List

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


def shorten_filenames(filenames: List[str]) -> List[str]:
    """
    Strip extensions, split on '_', remove common prefix/suffix tokens.

    Examples:
        ["Experiment_sample_A.mzML", "Experiment_sample_B.mzML"]
        -> ["sample_A", "sample_B"]
    """
    if not filenames:
        return []

    stems = [Path(f).stem for f in filenames]

    if len(stems) <= 1:
        return stems

    token_lists = [s.split("_") for s in stems]
    min_len = min(len(t) for t in token_lists)

    prefix_count = 0
    for i in range(min_len):
        if all(t[i] == token_lists[0][i] for t in token_lists):
            prefix_count += 1
        else:
            break

    suffix_count = 0
    for i in range(1, min_len - prefix_count + 1):
        if all(t[-i] == token_lists[0][-i] for t in token_lists):
            suffix_count += 1
        else:
            break

    result = []
    for tokens in token_lists:
        end = len(tokens) - suffix_count if suffix_count > 0 else len(tokens)
        remaining = tokens[prefix_count:end]
        if remaining:
            result.append("_".join(remaining))
        else:
            result.append("_".join(tokens))

    return result


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
    except (ValueError, FloatingPointError) as e:
        logger.warning(f"Entropy similarity failed: {e}")
        return float("nan")


def compute_shared_peak_dice(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Dice coefficient on shared annotated peaks."""
    _, _, shared = _get_matching_intensities(query_frame, reference_frame)
    ann_q = int(query_frame["is_annotated"].sum()) if "is_annotated" in query_frame.columns else 0
    ann_r = int(reference_frame["is_annotated"].sum()) if "is_annotated" in reference_frame.columns else 0
    denominator = ann_q + ann_r
    if denominator == 0:
        return float("nan")
    return float(2 * len(shared) / denominator)


def compute_shared_peak_jaccard(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Jaccard index on shared annotated peaks."""
    _, _, shared = _get_matching_intensities(query_frame, reference_frame)
    ann_q = int(query_frame["is_annotated"].sum()) if "is_annotated" in query_frame.columns else 0
    ann_r = int(reference_frame["is_annotated"].sum()) if "is_annotated" in reference_frame.columns else 0
    denominator = ann_q + ann_r - len(shared)
    if denominator == 0:
        return float("nan")
    return float(len(shared) / denominator)


def compute_shared_peak_ratio(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Shared peak ratio: shared / min(ann_q, ann_r)."""
    _, _, shared = _get_matching_intensities(query_frame, reference_frame)
    ann_q = int(query_frame["is_annotated"].sum()) if "is_annotated" in query_frame.columns else 0
    ann_r = int(reference_frame["is_annotated"].sum()) if "is_annotated" in reference_frame.columns else 0
    denominator = min(ann_q, ann_r)
    if denominator == 0:
        return float("nan")
    return float(len(shared) / denominator)


def compute_all_metrics(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> Dict[str, float]:
    """
    Compute all 7 comparison metrics plus context counts.

    Returns dict with keys:
        Metrics:  R2m, R2b, Cos, Ent, Dice, Jacc, ShR
        Counts:   total_q, total_r, ann_q, ann_r, shared
    """
    query_intensities, reference_intensities, shared_labels = _get_matching_intensities(query_frame, reference_frame)
    n_shared = len(shared_labels)

    query_annotated_count = int(query_frame["is_annotated"].sum()) if "is_annotated" in query_frame.columns else 0
    reference_annotated_count = int(reference_frame["is_annotated"].sum()) if "is_annotated" in reference_frame.columns else 0

    denom_dice = query_annotated_count + reference_annotated_count
    denom_jacc = query_annotated_count + reference_annotated_count - n_shared
    denom_shr = min(query_annotated_count, reference_annotated_count)

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
        "ann_q": query_annotated_count,
        "ann_r": reference_annotated_count,
        "shared": n_shared,
    }


def format_metrics_text(metrics: Dict[str, float]) -> str:
    """Format metrics as multi-line string for plot overlay."""
    def _fmt(val) -> str:
        if val is None or (isinstance(val, float) and math.isnan(val)):
            return "N/A"
        return f"{val:.3f}"

    return (
        f"R\u00b2 Match: {_fmt(metrics.get('R2m'))}\n"
        f"R\u00b2 Binned: {_fmt(metrics.get('R2b'))}\n"
        f"Cosine: {_fmt(metrics.get('Cos'))}\n"
        f"Entropy: {_fmt(metrics.get('Ent'))}\n"
        f"Dice: {_fmt(metrics.get('Dice'))}\n"
        f"Jaccard: {_fmt(metrics.get('Jacc'))}\n"
        f"Shared: {_fmt(metrics.get('ShR'))}"
    )


# Mapping from short keys to display column names (used by summary export)
METRIC_SHORT_TO_DISPLAY = {
    "R2m": "R\u00b2 match",
    "R2b": "R\u00b2 binned",
    "Cos": "Cosine",
    "Ent": "Entropy",
    "Dice": "Dice",
    "Jacc": "Jaccard",
    "ShR": "SharedRatio",
}
