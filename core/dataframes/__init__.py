"""Polars LazyFrames for openms_insight components.

Re-exports the public API of the submodules so callers can keep using
``from core.dataframes import ...`` without caring which file a name
lives in.
"""
from core.dataframes.psm_tables import (
    PSM_TABLE_SCHEMA,
    build_psm_table,
)
from core.dataframes.insight_frames import (
    HEATMAP_SCHEMA,
    PEAKS_SCHEMA,
    SEQUENCE_SCHEMA,
    build_heatmap_dataframe,
    build_peaks_dataframe,
    build_sequence_dataframe,
    build_visualization_frames,
)

__all__ = [
    "HEATMAP_SCHEMA",
    "PEAKS_SCHEMA",
    "PSM_TABLE_SCHEMA",
    "SEQUENCE_SCHEMA",
    "build_heatmap_dataframe",
    "build_peaks_dataframe",
    "build_psm_table",
    "build_sequence_dataframe",
    "build_visualization_frames",
]
