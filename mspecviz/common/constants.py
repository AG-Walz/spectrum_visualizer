"""Non-tunable constants used across the core package.

The tunable pipeline defaults (ion types, tolerances, Koina settings, export
options) live in :mod:`mspecviz.common.defaults`. This module is reserved for
constants that are not user-configurable — the session-state namespace and a
handful of plot-tuning knobs that ``mirror_plot`` / ``peak_annotation`` share.
"""

# Namespace for Streamlit session-state keys owned by the app. Used by
# ``streamlit_app.utils.state`` to avoid colliding with third-party widget keys.
STATE_PREFIX = "peptide_mirror"

# --- Plot tuning -----------------------------------------------------------
# Fraction of peak intensity used as vertical offset for annotation text.
ANN_OFFSET = 0.03

# Extra m/z padding on the plot x-axis (Da) so outer peaks aren't clipped.
MZ_WINDOW_PADDING = 50
