"""Shared helpers for Plotly + Matplotlib mirror plots.

Both backends start from an annotated DataFrame (``mz``, ``intensity``,
``annotation``) and need the same prep steps: add a default annotation column,
normalise intensities, and map ion types to colours. Centralising the logic
avoids colour-mapping drift between the two renderers.
"""
from typing import Dict

import pandas as pd

from mspecviz.common.ion_utils import ION_COLORS, get_ion_color
from mspecviz.common.normalization import normalize_intensities

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


def _prepare_spectrum_data(df: pd.DataFrame, normalize: bool = True,
                           color_map: Dict[str, str] = None) -> pd.DataFrame:
    """Prepare spectrum DataFrame for plotting (shared by Plotly and Matplotlib)."""
    df = _ensure_annotation_col(df.copy())
    if normalize:
        df = normalize_intensities(df)
    if color_map is None:
        color_map = ION_COLORS
    df["_color"] = df["annotation"].map(lambda a: get_ion_color(a, color_map))
    return df


def _get_annotated_mask(group: pd.DataFrame) -> pd.Series:
    """Return boolean mask for peaks that have a real annotation (not '' or 'unannotated')."""
    return group["annotation"].ne("unannotated") & group["annotation"].ne("")
