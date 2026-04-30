"""Shared helpers for Plotly + Matplotlib mirror plots.

Both backends start from an annotated DataFrame (``mz``, ``intensity``,
``annotation``) and need the same prep steps: add a default annotation column,
normalise intensities, and map ion types to colours. Centralising the logic
avoids colour-mapping drift between the two renderers.
"""
import pandas as pd

from mspecviz.common.fragment_utils import get_ion_color

# --- Plot tuning -----------------------------------------------------------
# Fraction of peak intensity used as vertical offset for annotation text.
ANN_OFFSET = 0.03

# Extra m/z padding on the plot x-axis (Da) so outer peaks aren't clipped.
MZ_WINDOW_PADDING = 50


def _ensure_annotation_col(df: pd.DataFrame) -> pd.DataFrame:
    """Return df with an 'annotation' column, adding 'unannotated' if absent."""
    if "annotation" not in df.columns:
        df = df.copy()
        df["annotation"] = "unannotated"
    return df


def _prepare_spectrum_data(df: pd.DataFrame, normalize: bool = True) -> pd.DataFrame:
    """Prepare spectrum DataFrame for plotting (shared by Plotly and Matplotlib)."""
    df = _ensure_annotation_col(df.copy())
    if normalize and not df.empty:
        max_intensity = df["intensity"].max()
        if max_intensity and max_intensity > 0:
            df["intensity"] = df["intensity"] / max_intensity
    df["_color"] = df["annotation"].map(get_ion_color)
    return df


def _get_annotated_mask(group: pd.DataFrame) -> pd.Series:
    """Return boolean mask for peaks that have a real annotation (not '' or 'unannotated')."""
    return group["annotation"].ne("unannotated") & group["annotation"].ne("")
