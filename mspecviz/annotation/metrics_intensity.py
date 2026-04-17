"""Intensity-based spectrum similarity metrics.

R\u00b2 (matched + binned) and cosine all rest on aligning intensities of
shared annotations. The shared helper :func:`_get_matching_intensities` lives
here; set-based metrics in :mod:`mspecviz.annotation.metrics_set` reuse it.
"""
import logging
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import linregress

logger = logging.getLogger(__name__)


def _get_matching_intensities(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Find overlapping annotated ions and return their intensity arrays.

    Returns:
        (query_intensities, reference_intensities, shared_labels)
    """
    query_annotated = query_frame[query_frame["is_annotated"]].copy()
    reference_annotated = reference_frame[reference_frame["is_annotated"]].copy()

    if query_annotated.empty or reference_annotated.empty:
        return np.array([]), np.array([]), []

    query_intensity_by_label = query_annotated.groupby("annotation")["intensity"].max().to_dict()
    reference_intensity_by_label = reference_annotated.groupby("annotation")["intensity"].max().to_dict()

    shared = sorted(set(query_intensity_by_label.keys()) & set(reference_intensity_by_label.keys()))
    if not shared:
        return np.array([]), np.array([]), []

    query_intensities = np.array([query_intensity_by_label[label] for label in shared])
    reference_intensities = np.array([reference_intensity_by_label[label] for label in shared])

    return query_intensities, reference_intensities, shared


def _r2_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """R^2 from pre-computed matched intensity arrays."""
    if len(query_intensities) < 2:
        return float("nan")
    _, _, r_value, _, _ = linregress(query_intensities, reference_intensities)
    return float(r_value ** 2)


def _cosine_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """Cosine similarity from pre-computed matched intensity arrays."""
    if len(query_intensities) == 0:
        return float("nan")
    dot = np.dot(query_intensities, reference_intensities)
    norm_query = np.linalg.norm(query_intensities)
    norm_reference = np.linalg.norm(reference_intensities)
    if norm_query == 0 or norm_reference == 0:
        return float("nan")
    return float(dot / (norm_query * norm_reference))


def compute_r2_matching(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """R^2 (matching) -- LinearRegression R^2 on matching annotated ion intensities."""
    query_intensities, reference_intensities, _ = _get_matching_intensities(query_frame, reference_frame)
    return _r2_from_arrays(query_intensities, reference_intensities)


def compute_r2_binned(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    n_bins: int = 1000,
) -> float:
    """R^2 (binned) -- LinearRegression R^2 on binned intensity profiles."""
    if query_frame.empty or reference_frame.empty:
        return float("nan")

    all_mz = np.concatenate([query_frame["mz"].values, reference_frame["mz"].values])
    mz_min = float(np.floor(all_mz.min()))
    mz_max = float(np.ceil(all_mz.max()))

    if mz_max <= mz_min:
        return float("nan")

    bin_edges = np.linspace(mz_min, mz_max, n_bins + 1)
    n_total = len(bin_edges) + 1
    query_bin_indices = np.digitize(query_frame["mz"].values, bin_edges)
    reference_bin_indices = np.digitize(reference_frame["mz"].values, bin_edges)
    query_sums = np.bincount(query_bin_indices, weights=query_frame["intensity"].values, minlength=n_total)
    reference_sums = np.bincount(reference_bin_indices, weights=reference_frame["intensity"].values, minlength=n_total)

    shared_mask = (query_sums > 0) & (reference_sums > 0)
    if shared_mask.sum() < 2:
        return float("nan")

    _, _, r_value, _, _ = linregress(query_sums[shared_mask], reference_sums[shared_mask])
    return float(r_value ** 2)


def compute_cosine_similarity(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Cosine similarity on matching annotated ion intensities."""
    query_intensities, reference_intensities, _ = _get_matching_intensities(query_frame, reference_frame)
    return _cosine_from_arrays(query_intensities, reference_intensities)
