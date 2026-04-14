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
from scipy.stats import linregress

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


def compute_spectral_entropy_similarity(
    query_df: pd.DataFrame,
    ref_df: pd.DataFrame,
) -> float:
    """Spectral entropy similarity (Li et al., Nature Methods 2021)."""
    if query_df.empty or ref_df.empty:
        return float("nan")

    if ms_entropy is None:
        logger.warning("ms_entropy not installed -- entropy similarity unavailable")
        return float("nan")

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
    """Dice coefficient on shared annotated peaks."""
    _, _, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0
    denominator = ann_q + ann_r
    if denominator == 0:
        return float("nan")
    return float(2 * len(shared) / denominator)


def compute_shared_peak_jaccard(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """Jaccard index on shared annotated peaks."""
    _, _, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0
    denominator = ann_q + ann_r - len(shared)
    if denominator == 0:
        return float("nan")
    return float(len(shared) / denominator)


def compute_shared_peak_ratio(query_df: pd.DataFrame, ref_df: pd.DataFrame) -> float:
    """Shared peak ratio: shared / min(ann_q, ann_r)."""
    _, _, shared = _get_matching_intensities(query_df, ref_df)
    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0
    denominator = min(ann_q, ann_r)
    if denominator == 0:
        return float("nan")
    return float(len(shared) / denominator)


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
    q_ints, r_ints, shared_labels = _get_matching_intensities(query_df, ref_df)
    n_shared = len(shared_labels)

    ann_q = int(query_df["is_annotated"].sum()) if "is_annotated" in query_df.columns else 0
    ann_r = int(ref_df["is_annotated"].sum()) if "is_annotated" in ref_df.columns else 0

    denom_dice = ann_q + ann_r
    denom_jacc = ann_q + ann_r - n_shared
    denom_shr = min(ann_q, ann_r)

    return {
        "R2m": _r2_from_arrays(q_ints, r_ints),
        "R2b": compute_r2_binned(query_df, ref_df),
        "Cos": _cosine_from_arrays(q_ints, r_ints),
        "Ent": compute_spectral_entropy_similarity(query_df, ref_df),
        "Dice": float(2 * n_shared / denom_dice) if denom_dice else float("nan"),
        "Jacc": float(n_shared / denom_jacc) if denom_jacc else float("nan"),
        "ShR": float(n_shared / denom_shr) if denom_shr else float("nan"),
        "total_q": len(query_df),
        "total_r": len(ref_df),
        "ann_q": ann_q,
        "ann_r": ann_r,
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
