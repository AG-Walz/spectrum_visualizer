"""Peak annotation + spectrum comparison metrics (re-exported flat for ``core.pipeline``)."""
from core.annotation.peak_matching import (
    annotate_spectrum,
    prepare_koina_spectrum_for_plot,
    prepare_spectrum_for_plot,
)
from core.annotation.intensity_metrics import (
    compute_cosine_similarity,
    compute_r2_binned,
    compute_r2_matching,
)
from core.annotation.comparison_metrics import (
    METRIC_SHORT_TO_DISPLAY,
    compute_all_metrics,
    compute_spectral_entropy_similarity,
    format_metrics_text,
)

__all__ = [
    # peak_matching
    "annotate_spectrum",
    "prepare_koina_spectrum_for_plot",
    "prepare_spectrum_for_plot",
    # intensity_metrics
    "compute_cosine_similarity",
    "compute_r2_binned",
    "compute_r2_matching",
    # comparison_metrics
    "METRIC_SHORT_TO_DISPLAY",
    "compute_all_metrics",
    "compute_spectral_entropy_similarity",
    "format_metrics_text",
]
