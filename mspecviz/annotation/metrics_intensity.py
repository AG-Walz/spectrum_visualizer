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
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Find overlapping annotated ions and return their intensity arrays.

    Returns:
        (query_intensities, ref_intensities, shared_labels)
    """
    q_ann = query_df[query_df["is_annotated"]].copy()
    r_ann = ref_df[ref_df["is_annotated"]].copy()

    if q_ann.empty or r_ann.empty:
        return np.array([]), np.array([]), []

    q_map = q_ann.groupby("annotation")["intensity"].max().to_dict()
    r_map = r_ann.groupby("annotation")["intensity"].max().to_dict()

    shared = sorted(set(q_map.keys()) & set(r_map.keys()))
    if not shared:
        return np.array([]), np.array([]), []

    q_ints = np.array([q_map[lbl] for lbl in shared])
    r_ints = np.array([r_map[lbl] for lbl in shared])

    return q_ints, r_ints, shared


def _r2_from_arrays(q_ints: np.ndarray, r_ints: np.ndarray) -> float:
    """R^2 from pre-computed matched intensity arrays."""
    if len(q_ints) < 2:
        return float("nan")
    _, _, r_value, _, _ = linregress(q_ints, r_ints)
    return float(r_value ** 2)


def _cosine_from_arrays(q_ints: np.ndarray, r_ints: np.ndarray) -> float:
    """Cosine similarity from pre-computed matched intensity arrays."""
    if len(q_ints) == 0:
        return float("nan")
    dot = np.dot(q_ints, r_ints)
    norm_q = np.linalg.norm(q_ints)
    norm_r = np.linalg.norm(r_ints)
    if norm_q == 0 or norm_r == 0:
        return float("nan")
    return float(dot / (norm_q * norm_r))


def compute_r2_matching(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """R^2 (matching) -- LinearRegression R^2 on matching annotated ion intensities."""
    q_ints, r_ints, _ = _get_matching_intensities(query_df, ref_df)
    return _r2_from_arrays(q_ints, r_ints)


def compute_r2_binned(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    n_bins: int = 1000,
) -> float:
    """R^2 (binned) -- LinearRegression R^2 on binned intensity profiles."""
    if query_df.empty or ref_df.empty:
        return float("nan")

    all_mz = np.concatenate([query_df["mz"].values, ref_df["mz"].values])
    mz_min = float(np.floor(all_mz.min()))
    mz_max = float(np.ceil(all_mz.max()))

    if mz_max <= mz_min:
        return float("nan")

    bin_edges = np.linspace(mz_min, mz_max, n_bins + 1)
    n_total = len(bin_edges) + 1
    q_idx = np.digitize(query_df["mz"].values, bin_edges)
    r_idx = np.digitize(ref_df["mz"].values, bin_edges)
    q_sums = np.bincount(q_idx, weights=query_df["intensity"].values, minlength=n_total)
    r_sums = np.bincount(r_idx, weights=ref_df["intensity"].values, minlength=n_total)

    shared_mask = (q_sums > 0) & (r_sums > 0)
    if shared_mask.sum() < 2:
        return float("nan")

    _, _, r_value, _, _ = linregress(q_sums[shared_mask], r_sums[shared_mask])
    return float(r_value ** 2)


def compute_cosine_similarity(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """Cosine similarity on matching annotated ion intensities."""
    q_ints, r_ints, _ = _get_matching_intensities(query_df, ref_df)
    return _cosine_from_arrays(q_ints, r_ints)
