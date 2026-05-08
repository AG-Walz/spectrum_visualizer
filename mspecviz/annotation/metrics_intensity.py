"""Intensity-based spectrum similarity metrics.

R² (matched + binned) and cosine all rest on aligning intensities of
shared annotations. The shared helper :func:`_get_matching_intensities` lives
here; set-based metrics in :mod:`mspecviz.annotation.comparison_metrics` reuse it.
"""
import logging

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine as _cosine_distance
from scipy.stats import linregress

logger = logging.getLogger(__name__)


def _get_matching_intensities(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Find overlapping annotated ions; returns ``(q_int, r_int, shared_labels)``.

    Within each spectrum, multiple peaks may share an annotation (e.g. isotopes
    or charge states with the same string label). We collapse to the brightest
    per label, then keep only labels present in both spectra.
    """
    query_max = query_frame[query_frame["is_annotated"]].groupby("annotation")["intensity"].max()
    reference_max = reference_frame[reference_frame["is_annotated"]].groupby("annotation")["intensity"].max()
    shared = sorted(set(query_max.index) & set(reference_max.index))
    if not shared:
        return np.array([]), np.array([]), []
    return query_max.reindex(shared).to_numpy(), reference_max.reindex(shared).to_numpy(), shared


def _r2_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """R^2 from pre-computed matched intensity arrays."""
    if len(query_intensities) < 2:
        return float("nan")
    return float(linregress(query_intensities, reference_intensities).rvalue ** 2)


def _cosine_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """Cosine similarity from pre-computed matched intensity arrays."""
    if len(query_intensities) == 0 or not np.any(query_intensities) or not np.any(reference_intensities):
        return float("nan")
    return 1.0 - float(_cosine_distance(query_intensities, reference_intensities))


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

    mz_min = float(np.floor(min(query_frame["mz"].min(), reference_frame["mz"].min())))
    mz_max = float(np.ceil(max(query_frame["mz"].max(), reference_frame["mz"].max())))
    if mz_max <= mz_min:
        return float("nan")

    bin_edges = np.linspace(mz_min, mz_max, n_bins + 1)
    n_total = len(bin_edges) + 1
    query_bins = np.digitize(query_frame["mz"].to_numpy(), bin_edges)
    reference_bins = np.digitize(reference_frame["mz"].to_numpy(), bin_edges)
    query_sums = np.bincount(query_bins, weights=query_frame["intensity"].to_numpy(), minlength=n_total)
    reference_sums = np.bincount(reference_bins, weights=reference_frame["intensity"].to_numpy(), minlength=n_total)

    shared_mask = (query_sums > 0) & (reference_sums > 0)
    if shared_mask.sum() < 2:
        return float("nan")
    return float(linregress(query_sums[shared_mask], reference_sums[shared_mask]).rvalue ** 2)


def compute_cosine_similarity(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Cosine similarity on matching annotated ion intensities."""
    query_intensities, reference_intensities, _ = _get_matching_intensities(query_frame, reference_frame)
    return _cosine_from_arrays(query_intensities, reference_intensities)
