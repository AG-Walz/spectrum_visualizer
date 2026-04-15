"""
Intensity normalization helpers.

Both ``mspecviz.plotting`` and ``mspecviz.annotation`` need the same
"divide intensity by its max, unless max is zero" pattern. Centralising the
check means a fix to the degenerate-input case (empty df / all-zero) lands in
one place rather than drifting between modules.
"""
from __future__ import annotations

import pandas as pd


def normalize_intensities(df: pd.DataFrame, column: str = "intensity") -> pd.DataFrame:
    """Divide ``df[column]`` by its max in place and return the same df.

    Safe when ``df`` is empty or the column is all-zero: the column is left
    unchanged. The caller is responsible for passing a df they own (we mutate).
    """
    if df.empty:
        return df
    max_val = df[column].max()
    if max_val and max_val > 0:
        df[column] = df[column] / max_val
    return df
