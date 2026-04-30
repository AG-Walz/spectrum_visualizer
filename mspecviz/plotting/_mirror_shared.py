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

# Colour palette keyed by ``ion_type`` bucket (see ``annotation.peak_matching.get_ion_type``).
ION_COLORS = {
    "a": "#FFCC00",
    "b": "#0066CC",
    "c": "#00CC66",
    "x": "#CC6600",
    "y": "#CC0000",
    "z": "#6600CC",
    "precursor": "#8B008B",
    "immonium": "#00CED1",
    "unannotated": "#AAAAAA",
}


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
