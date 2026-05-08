"""Peak annotation and spectrum comparison metrics.

Re-exports the public API of the three submodules so callers can keep using
``from mspecviz.annotation import ...`` without caring which file a name lives in.
"""
from mspecviz.annotation.peak_matching import (
    annotate_spectrum,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from mspecviz.annotation.metrics_intensity import (
    compute_cosine_similarity,
    compute_r2_binned,
    compute_r2_matching,
)
from mspecviz.annotation.comparison_metrics import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    compute_shared_peak_dice,
    compute_shared_peak_jaccard,
    compute_shared_peak_ratio,
    compute_spectral_entropy_similarity,
    format_metrics_text,
    shorten_filenames,
)

__all__ = [
    # peak_matching
    "annotate_spectrum",
    "prepare_koina_spectrum_for_plot",
    "prepare_spectrum_for_plot",
    # metrics_intensity
    "compute_cosine_similarity",
    "compute_r2_binned",
    "compute_r2_matching",
    # comparison_metrics
    "METRIC_SHORT_TO_DISPLAY",
    "compute_all_metrics",
    "compute_shared_peak_dice",
    "compute_shared_peak_jaccard",
    "compute_shared_peak_ratio",
    "compute_spectral_entropy_similarity",
    "format_metrics_text",
    "shorten_filenames",
]
