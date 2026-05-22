"""Compute 3 intensity-based similarity scores (cosine, R² matched, R² binned) 
between two annotated spectra (used by comparison_metrics.compute_all_metrics).
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
) -> tuple[np.ndarray, np.ndarray]:
    """Aligns 2 spectra on the ions they both annotated returning their intensities as parallel arrays. 
    Each array is ordered by annotation and has length equal to the number of shared annotations. 
    Returns a pair of arrays in the form (query, reference) of length n-shared.
    """
    # Keep only annotated peaks, groups by the annotation label, and takes .max() intensity per label. 
    # This collapses duplicate labels — isotopes or alternative charge states that share a label collapse to their brightest peak.
    query_max = query_frame[query_frame["is_annotated"]].groupby("annotation")["intensity"].max()
    reference_max = reference_frame[reference_frame["is_annotated"]].groupby("annotation")["intensity"].max()
    # Intersect the two label sets and reindex both to the same sorted shared labels.
    shared = sorted(set(query_max.index) & set(reference_max.index))
    # Empty arrays if nothing is shared.
    if not shared:
        return np.array([]), np.array([])
    # Returns two equal-length, position-aligned NumPy arrays (query, reference). 
    return query_max.reindex(shared).to_numpy(), reference_max.reindex(shared).to_numpy()


def _r2_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """Squared Pearson correlation (R²) from already-aligned intensity arrays; NaN with fewer than 2 points.
    """
    if len(query_intensities) < 2:
        return float("nan")
    return float(linregress(query_intensities, reference_intensities).rvalue ** 2)


def _cosine_from_arrays(query_intensities: np.ndarray, reference_intensities: np.ndarray) -> float:
    """Cosine similarity from already-aligned intensity arrays; NaN if either side is all-zero.
    """
    if len(query_intensities) == 0 or not np.any(query_intensities) or not np.any(reference_intensities):
        return float("nan")
    return 1.0 - float(_cosine_distance(query_intensities, reference_intensities))


def compute_r2_matching(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """R² on intensities of ions annotated in both spectra (R2m column).
    """
    query_intensities, reference_intensities = _get_matching_intensities(query_frame, reference_frame)
    return _r2_from_arrays(query_intensities, reference_intensities)


def compute_r2_binned(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    n_bins: int = 1000,
) -> float:
    """R² on intensity profiles after binning both spectra onto a shared m/z grid (R2b column).
    """ 
    # Avoid binning empty inputs or a zero-width m/z range 
    if query_frame.empty or reference_frame.empty:
        return float("nan")

    # Define a shared m/z grid that covers both spectra, with a fixed number of bins. The grid edges are rounded to the nearest integer.
    mz_min = float(np.floor(min(query_frame["mz"].min(), reference_frame["mz"].min())))
    mz_max = float(np.ceil(max(query_frame["mz"].max(), reference_frame["mz"].max())))
    if mz_max <= mz_min:
        return float("nan")

    # Bin spectra: sum intensities into a shared 1000-bin grid spanning the combined m/z range
    # using all bins, not just those with annotations, to capture overall profile similarity.
    bin_edges = np.linspace(mz_min, mz_max, n_bins + 1)
    query_sums, _ = np.histogram(query_frame["mz"], bins=bin_edges, weights=query_frame["intensity"])
    reference_sums, _ = np.histogram(reference_frame["mz"], bins=bin_edges, weights=reference_frame["intensity"])

    # Regress: only bins with signal in both spectra count toward the fit; NaN if fewer than 2 such bins.
    shared_mask = (query_sums > 0) & (reference_sums > 0)
    if shared_mask.sum() < 2:
        return float("nan")
    return float(linregress(query_sums[shared_mask], reference_sums[shared_mask]).rvalue ** 2)


def compute_cosine_similarity(query_frame: pd.DataFrame, reference_frame: pd.DataFrame) -> float:
    """Cosine similarity on intensities of ions annotated in both spectra (Cos column).
    """
    query_intensities, reference_intensities = _get_matching_intensities(query_frame, reference_frame)
    return _cosine_from_arrays(query_intensities, reference_intensities)
