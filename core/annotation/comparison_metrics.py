"""Spectrum pair scoring layer.

Given two already-annotated spectra (a query and a reference of the same
peptide), produces:

- Three intensity-based similarity scores (cosine, R² matched, R² binned),
  computed from the shared-annotation intensity vectors.
- A spectral-entropy similarity (Li et al., Nature Methods 2021).
- `compute_all_metrics`, which combines those with three set-overlap metrics
  (Dice, Jaccard, SharedRatio) into the seven-metric row used by the
  comparison table, export spreadsheet, and mirror-plot corner overlay.

Called by core.pipeline.query_processing and core.plotting.plot_worker.
"""
import logging
import math

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine as _cosine_distance
from scipy.stats import linregress

from core.defaults import DEFAULTS

try:
    import ms_entropy
except ImportError:
    ms_entropy = None

logger = logging.getLogger(__name__)


# ============================================================================
# Intensity-aligned helpers (private — used by the R²/cosine scorers below)
# ============================================================================

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


# ============================================================================
# Public intensity scorers
# ============================================================================

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


# ============================================================================
# Spectral-entropy similarity + composite metric row
# ============================================================================

def compute_spectral_entropy_similarity(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    tolerance_da: float = DEFAULTS.mz_tolerance,
) -> float:
    """Wraper for spectral entropy similarity between two annotated spectra
    (Li et al., Nature Methods 2021)
    """
    # Return NaN for empty spectra
    if query_frame.empty or reference_frame.empty:
        return float("nan")
    # Return NaN if ms_entropy is not installed
    if ms_entropy is None:
        logger.warning("ms_entropy not installed -- entropy similarity unavailable")
        return float("nan")

    # Pack spectrum in 2 columns (mz, intensity) array
    query_spectrum_array = np.column_stack([query_frame["mz"].values, query_frame["intensity"].values])
    reference_spectrum_array = np.column_stack([reference_frame["mz"].values, reference_frame["intensity"].values])

    # Compute: downgrade numeric failures to NaN so the comparison table never breaks
    try:
        similarity = ms_entropy.calculate_entropy_similarity(
            query_spectrum_array, reference_spectrum_array, ms2_tolerance_in_da=tolerance_da
        )
        return float(similarity)
    except (ValueError, FloatingPointError) as error:
        logger.warning(f"Entropy similarity failed: {error}")
        return float("nan")


def compute_all_metrics(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    tolerance_da: float = DEFAULTS.mz_tolerance,
) -> dict[str, float | int]:
    """Compute one row of the comparison table: 7 similarity metrics (3 Intensity based: R2m, R2b, Cos,
    1 Entropy-based: Ent, 3 Ion-matching: Dice, Jacc, ShR) plus peak/annotation counts.
    """
    # Align the query and reference spectra on their shared annotated ions
    query_intensities, reference_intensities = _get_matching_intensities(query_frame, reference_frame)
    n_shared = len(query_intensities)

    # Count unique ion labels per spectrum (duplicate isotopes/charges collapse to one)
    query_annotation_count = int(query_frame[query_frame["is_annotated"]]["annotation"].nunique())
    reference_annotation_count = int(reference_frame[reference_frame["is_annotated"]]["annotation"].nunique())

    # Denominators for the three set-overlap metrics (Dice, Jaccard, SharedRatio)
    denom_dice = query_annotation_count + reference_annotation_count
    denom_jacc = query_annotation_count + reference_annotation_count - n_shared
    denom_shr = min(query_annotation_count, reference_annotation_count)

    # Assemble: 7 metrics + the context counts the export and UI table show alongside them
    return {
        "R2m": _r2_from_arrays(query_intensities, reference_intensities),
        "R2b": compute_r2_binned(query_frame, reference_frame),
        "Cos": _cosine_from_arrays(query_intensities, reference_intensities),
        "Ent": compute_spectral_entropy_similarity(query_frame, reference_frame, tolerance_da=tolerance_da),
        "Dice": float(2 * n_shared / denom_dice) if denom_dice else float("nan"),
        "Jacc": float(n_shared / denom_jacc) if denom_jacc else float("nan"),
        "ShR": float(n_shared / denom_shr) if denom_shr else float("nan"),
        "query_peak_count": len(query_frame),
        "reference_peak_count": len(reference_frame),
        "query_annotation_count": query_annotation_count,
        "reference_annotation_count": reference_annotation_count,
        "shared": n_shared,
    }


def format_metrics_text(metrics: dict[str, float | int]) -> str:
    """Format the 7 metrics as a multi-line string for the mirror-plot corner overlay.
    """
    # Format None and NaN as "N/A"
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


# Short metric keys -> human-readable column names used in the export spreadsheet
METRIC_SHORT_TO_DISPLAY = {
    "R2m": "R² match",
    "R2b": "R² binned",
    "Cos": "Cosine",
    "Ent": "Entropy",
    "Dice": "Dice",
    "Jacc": "Jaccard",
    "ShR": "SharedRatio",
}