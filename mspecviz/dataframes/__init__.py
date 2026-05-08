"""Polars LazyFrames for openms_insight components.

Re-exports the public API of the submodules so callers can keep using
``from mspecviz.dataframes import ...`` without caring which file a name
lives in.
"""
from mspecviz.dataframes.psm_tables import (
    PSM_TABLE_SCHEMA,
    _enforce_schema,
    build_psm_table,
    validate_dataframe,
)
from mspecviz.dataframes.insight_frames import (
    HEATMAP_SCHEMA,
    PEAKS_SCHEMA,
    SEQUENCE_SCHEMA,
    _proforma_to_openms_str,
    build_heatmap_dataframe,
    build_peaks_dataframe,
    build_sequence_dataframe,
    save_to_parquet,
)

__all__ = [
    # psm_tables
    "PSM_TABLE_SCHEMA",
    "_enforce_schema",
    "build_psm_table",
    "validate_dataframe",
    # insight_frames
    "HEATMAP_SCHEMA",
    "PEAKS_SCHEMA",
    "SEQUENCE_SCHEMA",
    "_proforma_to_openms_str",
    "build_heatmap_dataframe",
    "build_peaks_dataframe",
    "build_sequence_dataframe",
    "save_to_parquet",
]
