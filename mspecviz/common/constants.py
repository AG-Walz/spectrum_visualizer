"""
Centralised fallback defaults + tunable constants.

Project config lives in JSON (``streamlit_app/config/settings.json``,
``streamlit_app/config/default-parameters.json``) and is the primary source of
truth. This module exposes a single ``DEFAULTS`` dict with the last-resort
fallback values used when neither JSON file — nor an active workspace
``params.json`` — supplies the key. Keep the three in sync.

Plot tuning constants live here too so ``mirror_plot`` / ``peak_annotation``
don't drift on shared visual choices.
"""

# --- Workspace / session-state ---------------------------------------------
# Namespace for session-state keys owned by the Streamlit app. Used by
# ``streamlit_app.utils.state`` to avoid colliding with third-party widget keys.
STATE_PREFIX = "peptide_mirror"

# --- Pipeline defaults -----------------------------------------------------
DEFAULTS: dict = {
    # Selection
    "selected-idxmls": [],
    "selected-mzmls": [],
    "target-peptides": "",
    "reference-substring": "",

    # Annotation
    "ion-types": ["b", "y"],
    "mz-tolerance": 0.02,            # Da
    "add-losses": False,
    "add-first-prefix-ion": True,
    "remove-precursor-peaks": False,
    "intensity-mode": "Relative",

    # Performance
    "workers": 4,

    # Koina AI prediction
    "koina-enabled": False,
    "koina-model": "ms2pip_timsTOF2024",
    "koina-ce": 35,
    "koina-instrument": "TIMSTOF",

    # UI / export
    "image-format": "png",
    "export-dir": "",
}

# --- Plot tuning -----------------------------------------------------------
# Fraction of peak intensity used as vertical offset for annotation text.
ANN_OFFSET = 0.03

# Extra m/z padding on the plot x-axis (Da) so outer peaks aren't clipped.
MZ_WINDOW_PADDING = 50
