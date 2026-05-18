"""Full 7-metric comparison row + entropy + plot-overlay text (called by ``core.pipeline.query_processing``; 17 ion-stats fields added later by ``core.plotting.plot_worker``)."""

import logging
import math
import numpy as np
import pandas as pd

from core.annotation.metrics_intensity import (
    _cosine_from_arrays,
    _get_matching_intensities,
    _r2_from_arrays,
    compute_r2_binned,
)
from core.common.defaults import DEFAULTS

try:
    import ms_entropy
except ImportError:
    ms_entropy = None

logger = logging.getLogger(__name__)


def compute_spectral_entropy_similarity(
    query_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    tolerance_da: float = DEFAULTS.mz_tolerance,
) -> float:
    """Spectral entropy similarity between two annotated spectra (Li et al., Nature Methods 2021); NaN if empty or ``ms_entropy`` is missing."""
    # Guard: empty spectra or optional dependency unavailable -> NaN so the row still tabulates
    if query_frame.empty or reference_frame.empty:
        return float("nan")

    if ms_entropy is None:
        logger.warning("ms_entropy not installed -- entropy similarity unavailable")
        return float("nan")

    # Pack: ms_entropy expects each spectrum as a 2-column (mz, intensity) array
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
    """One row of the comparison table: 7 similarity metrics (R2m, R2b, Cos, Ent, Dice, Jacc, ShR) plus peak/annotation counts."""
    # Align: shared annotated ions feed the three intensity-based metrics
    query_intensities, reference_intensities = _get_matching_intensities(query_frame, reference_frame)
    n_shared = len(query_intensities)

    # Count: unique ion labels per spectrum (duplicate isotopes/charges collapse to one)
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
    """Format the 7 metrics as a multi-line string for the mirror-plot corner overlay."""
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
