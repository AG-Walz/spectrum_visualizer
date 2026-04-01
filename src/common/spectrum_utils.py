"""
Shared spectrum utilities for MS/MS data handling.

Used by src.mirror_plot and src.peak_annotation.
"""
import pandas as pd


def normalize_intensities(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize df['intensity'] to [0, 1] in-place and return df."""
    max_int = df["intensity"].max()
    if max_int > 0:
        df["intensity"] = df["intensity"] / max_int
    return df


def ensure_annotation_col(df: pd.DataFrame) -> pd.DataFrame:
    """Return df with an 'annotation' column, adding 'unannotated' if absent."""
    if "annotation" not in df.columns:
        df = df.copy()
        df["annotation"] = "unannotated"
    return df
