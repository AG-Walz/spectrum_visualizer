"""
Spectrum comparison metrics for MS/MS spectra.

Pure Python backend module (no Streamlit). Computes quantitative similarity
metrics between query and reference annotated spectra.
"""
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Filename shortening
# ---------------------------------------------------------------------------

def shorten_filenames(filenames: List[str]) -> List[str]:
    """
    Strip extensions, split on '_', remove common prefix/suffix tokens.

    Examples:
        ["ProCa_Anna_sample_A.mzML", "ProCa_Anna_sample_B.mzML"]
        -> ["sample_A", "sample_B"]

        Single file -> stem (no shortening needed)
        Empty list  -> empty list
    """
    if not filenames:
        return []

    # Strip extensions to get stems
    stems = [Path(f).stem for f in filenames]

    if len(stems) <= 1:
        return stems

    # Split each stem into tokens on '_'
    token_lists = [s.split("_") for s in stems]

    # Find common prefix tokens
    min_len = min(len(t) for t in token_lists)
    prefix_count = 0
    for i in range(min_len):
        if all(t[i] == token_lists[0][i] for t in token_lists):
            prefix_count += 1
        else:
            break

    # Find common suffix tokens
    suffix_count = 0
    for i in range(1, min_len - prefix_count + 1):
        if all(t[-i] == token_lists[0][-i] for t in token_lists):
            suffix_count += 1
        else:
            break

    # Strip common prefix and suffix
    result = []
    for tokens in token_lists:
        end = len(tokens) - suffix_count if suffix_count > 0 else len(tokens)
        remaining = tokens[prefix_count:end]
        if remaining:
            result.append("_".join(remaining))
        else:
            # All tokens were common — fall back to original stem
            result.append("_".join(tokens))

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_matching_intensities(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Find overlapping annotated ions and return their intensity arrays.

    Filters both DataFrames to is_annotated == True, finds common annotation
    labels, returns ordered intensity arrays for shared annotations.

    Returns:
        (query_intensities, ref_intensities, shared_labels)
        Arrays may be empty if no overlap.
    """
    q_ann = query_df[query_df["is_annotated"] == True].copy()
    r_ann = ref_df[ref_df["is_annotated"] == True].copy()

    if q_ann.empty or r_ann.empty:
        return np.array([]), np.array([]), []

    # Build label -> intensity maps (take max intensity if duplicate annotations)
    q_map: Dict[str, float] = {}
    for _, row in q_ann.iterrows():
        label = row["annotation"]
        if label not in q_map or row["intensity"] > q_map[label]:
            q_map[label] = row["intensity"]

    r_map: Dict[str, float] = {}
    for _, row in r_ann.iterrows():
        label = row["annotation"]
        if label not in r_map or row["intensity"] > r_map[label]:
            r_map[label] = row["intensity"]

    shared = sorted(set(q_map.keys()) & set(r_map.keys()))
    if not shared:
        return np.array([]), np.array([]), []

    q_ints = np.array([q_map[lbl] for lbl in shared])
    r_ints = np.array([r_map[lbl] for lbl in shared])

    return q_ints, r_ints, shared


# ---------------------------------------------------------------------------
# Metric functions — all take (query_df, ref_df) annotated DataFrames
# ---------------------------------------------------------------------------

def compute_r2_matching(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """
    R^2 (matching) — LinearRegression R^2 on matching annotated ion intensities.

    Reference: specviz spectrum_correlation()
    """
    q_ints, r_ints, shared = _get_matching_intensities(query_df, ref_df)
    if len(shared) < 2:
        return float("nan")

    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(q_ints.reshape(-1, 1), r_ints)
    return float(model.score(q_ints.reshape(-1, 1), r_ints))


def compute_r2_binned(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    n_bins: int = 1000,
) -> float:
    """
    R^2 (binned) — LinearRegression R^2 on binned intensity profiles.

    Bins the full m/z range into `n_bins` bins and sums intensities per bin
    for both spectra, then fits a linear regression on the *intersection*
    of bins (where both spectra have signal).

    Reference: specviz new_spectrum_correlation()
    """
    if query_df.empty or ref_df.empty:
        return float("nan")

    all_mz = np.concatenate([query_df["mz"].values, ref_df["mz"].values])
    mz_min = float(np.floor(all_mz.min()))
    mz_max = float(np.ceil(all_mz.max()))

    if mz_max <= mz_min:
        return float("nan")

    bin_edges = np.linspace(mz_min, mz_max, n_bins + 1)

    # Use np.digitize (matching specviz approach) and sum per bin
    from collections import defaultdict

    def _bin_intensities(mz_arr, int_arr):
        bin_dict = defaultdict(float)
        indices = np.digitize(mz_arr, bin_edges)
        for b, intensity in zip(indices, int_arr):
            bin_dict[b] += intensity
        return bin_dict

    q_bins = _bin_intensities(query_df["mz"].values, query_df["intensity"].values)
    r_bins = _bin_intensities(ref_df["mz"].values, ref_df["intensity"].values)

    # Intersection of bins (where both spectra have signal) — specviz convention
    shared_bins = sorted(set(q_bins.keys()) & set(r_bins.keys()))
    if len(shared_bins) < 2:
        return float("nan")

    q_arr = np.array([q_bins[b] for b in shared_bins])
    r_arr = np.array([r_bins[b] for b in shared_bins])

    if np.all(q_arr == 0) or np.all(r_arr == 0):
        return float("nan")

    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(q_arr.reshape(-1, 1), r_arr)
    return float(model.score(q_arr.reshape(-1, 1), r_arr))


def compute_cosine_similarity(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """
    Cosine similarity on matching annotated ion intensities.

    cos = dot(q, r) / (||q|| * ||r||)

    Reference: Stein & Scott 1994
    """
    q_ints, r_ints, shared = _get_matching_intensities(query_df, ref_df)
    if len(shared) == 0:
        return float("nan")

    dot = np.dot(q_ints, r_ints)
    norm_q = np.linalg.norm(q_ints)
    norm_r = np.linalg.norm(r_ints)

    if norm_q == 0 or norm_r == 0:
        return float("nan")

    return float(dot / (norm_q * norm_r))


def compute_spectral_entropy_similarity(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
) -> float:
    """
    Spectral entropy similarity using the ms_entropy library.

    S = 1 - (2*H(merged) - H(q) - H(r)) / ln(4)

    Reference: Li et al., Nature Methods 2021
    """
    if query_df.empty or ref_df.empty:
        return float("nan")

    try:
        import ms_entropy
    except ImportError:
        logger.warning("ms_entropy not installed — entropy similarity unavailable")
        return float("nan")

    # ms_entropy expects 2D arrays of shape (n_peaks, 2) with [mz, intensity]
    q_spec = np.column_stack([query_df["mz"].values, query_df["intensity"].values])
    r_spec = np.column_stack([ref_df["mz"].values, ref_df["intensity"].values])

    try:
        similarity = ms_entropy.calculate_entropy_similarity(
            q_spec, r_spec, ms2_tolerance_in_da=0.02
        )
        return float(similarity)
    except Exception as e:
        logger.warning(f"Entropy similarity failed: {e}")
        return float("nan")


def compute_shared_peak_dice(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """
    Dice coefficient on shared annotated peaks.

    Dice = 2 * shared / (ann_q + ann_r)
    """
    q_ints, r_ints, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = query_df["is_annotated"].sum() if "is_annotated" in query_df.columns else 0
    ann_r = ref_df["is_annotated"].sum() if "is_annotated" in ref_df.columns else 0

    denominator = ann_q + ann_r
    if denominator == 0:
        return float("nan")

    return float(2 * len(shared) / denominator)


def compute_shared_peak_jaccard(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """
    Jaccard index on shared annotated peaks.

    Jaccard = shared / (ann_q + ann_r - shared)
    """
    q_ints, r_ints, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0

    denominator = ann_q + ann_r - len(shared)
    if denominator == 0:
        return float("nan")

    return float(len(shared) / denominator)


def compute_shared_peak_ratio(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """
    Shared peak ratio: shared / min(ann_q, ann_r).

    Measures completeness relative to the smaller set.
    """
    q_ints, r_ints, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0

    denominator = min(ann_q, ann_r)
    if denominator == 0:
        return float("nan")

    return float(len(shared) / denominator)


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def compute_all_metrics(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
) -> Dict[str, float]:
    """
    Compute all 7 comparison metrics plus context counts.

    Returns dict with keys:
        Metrics:  R2m, R2b, Cos, Ent, Dice, Jacc, ShR
        Counts:   total_q, total_r, ann_q, ann_r, shared
    """
    # Pre-compute shared info (avoids redundant work across metrics)
    _, _, shared_labels = _get_matching_intensities(query_df, ref_df)

    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0

    return {
        "R2m": compute_r2_matching(query_df, ref_df),
        "R2b": compute_r2_binned(query_df, ref_df),
        "Cos": compute_cosine_similarity(query_df, ref_df),
        "Ent": compute_spectral_entropy_similarity(query_df, ref_df),
        "Dice": compute_shared_peak_dice(query_df, ref_df),
        "Jacc": compute_shared_peak_jaccard(query_df, ref_df),
        "ShR": compute_shared_peak_ratio(query_df, ref_df),
        "total_q": len(query_df),
        "total_r": len(ref_df),
        "ann_q": ann_q,
        "ann_r": ann_r,
        "shared": len(shared_labels),
    }


def format_metrics_text(metrics: Dict[str, float]) -> str:
    """
    Format metrics as compact multi-line string for plot overlay.

    Shows top 5 metrics (R2m, R2b, Cos, Ent, Dice) in a 2-2-1 layout.
    NaN values shown as "N/A".
    """
    def _fmt(key: str) -> str:
        val = metrics.get(key, float("nan"))
        if val is None or (isinstance(val, float) and math.isnan(val)):
            return f"{key}=N/A"
        return f"{key}={val:.2f}"

    line1 = f"{_fmt('R2m')}  {_fmt('R2b')}"
    line2 = f"{_fmt('Cos')}  {_fmt('Ent')}"
    line3 = _fmt("Dice")

    return f"{line1}\n{line2}\n{line3}"
