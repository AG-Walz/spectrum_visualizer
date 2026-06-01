"""Peak annotation + spectrum comparison metrics (re-exported flat for core.pipeline).
"""
from core.annotation.peak_matching import (
    annotate_spectrum,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.annotation.comparison_metrics import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    compute_r2_binned,
    compute_spectral_entropy_similarity,
    format_metrics_text,
)

__all__ = [
    # peak_matching
    "annotate_spectrum",
    "prepare_koina_spectrum_for_plot",
    "prepare_spectrum_for_plot",
    # comparison_metrics
    "METRIC_SHORT_TO_DISPLAY",
    "compute_all_metrics",
    "compute_r2_binned",
    "compute_spectral_entropy_similarity",
    "format_metrics_text",
]