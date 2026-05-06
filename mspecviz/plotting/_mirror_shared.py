"""Shared helpers for Plotly + Matplotlib mirror plots.

Both backends start from an annotated DataFrame (``mz``, ``intensity``,
``annotation``, ``ion_type``) and need the same prep steps: add default
columns when absent, normalise intensities, map ion types to colours.
Centralising avoids colour-mapping drift between the two renderers.
"""
import pandas as pd

# --- Plot tuning -----------------------------------------------------------
# Fraction of peak intensity used as vertical offset for annotation text.
ANN_OFFSET = 0.03

# Extra m/z padding on the plot x-axis (Da) so outer peaks aren't clipped.
MZ_WINDOW_PADDING = 50

# Colour palette keyed by ``ion_type`` bucket (see
# ``annotation.peak_matching.get_ion_type``).
#
# These are the ``spectrum_utils.plot.colors`` defaults (Material Design
# 700-tone) — the de-facto standard in the Python MS plotting ecosystem and
# the same palette the upstream annotator (``MsmsSpectrum.annotate_proforma``)
# uses internally. Cross-checked against MaxQuant Viewer and OpenMS
# TOPPView, which agree on the b=blue / y=red / a=green / immonium=grey
# choices.
ION_COLORS = {
    "a":           "#388E3C",  # green       (Material Green 700)
    "b":           "#1976D2",  # blue        (Material Blue 700)
    "c":           "#00796B",  # teal        (Material Teal 700)
    "x":           "#7B1FA2",  # purple      (Material Purple 700)
    "y":           "#D32F2F",  # red         (Material Red 700)
    "z":           "#F57C00",  # orange      (Material Orange 700)
    "internal":    "#FBC02D",  # yellow      (Material Yellow 700) — spectrum_utils' "m"
    "precursor":   "#512DA8",  # deep purple (Material Deep Purple 700)
    "immonium":    "#455A64",  # blue-grey   (Material Blue Grey 700)
    "unannotated": "#9E9E9E",  # grey        (Material Grey 500)
}
# Note: neutral-loss peaks (``b3-H2O``, ``y5-NH3``) intentionally don't get a
# bucket of their own — they keep the parent ion's colour, matching the
# spectrum_utils plotting convention.


def _ensure_plot_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return df with the columns the prep pipeline expects, adding ``"unannotated"`` fallbacks."""
    for col in ("annotation", "ion_type"):
        if col not in df.columns:
            df = df.copy()
            df[col] = "unannotated"
    return df


def _prepare_spectrum_data(df: pd.DataFrame, normalize: bool = True) -> pd.DataFrame:
    """Prepare spectrum DataFrame for plotting (shared by Plotly and Matplotlib)."""
    df = _ensure_plot_columns(df.copy())
    if normalize and not df.empty:
        max_intensity = df["intensity"].max()
        if max_intensity and max_intensity > 0:
            df["intensity"] = df["intensity"] / max_intensity
    df["_color"] = df["ion_type"].map(ION_COLORS).fillna(ION_COLORS["unannotated"])
    return df


def _get_annotated_mask(group: pd.DataFrame) -> pd.Series:
    """Return boolean mask for peaks that have a real annotation (not '' or 'unannotated')."""
    return group["annotation"].ne("unannotated") & group["annotation"].ne("")
