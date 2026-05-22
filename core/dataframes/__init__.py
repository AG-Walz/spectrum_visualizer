"""Polars LazyFrames feeding openms_insight Table/LinePlot/SequenceView/Heatmap 
(imported by cli.main and streamlit_app.views.load_annotate.load_annotate_actions).
"""
# The PSM selector Table — the row-per-PSM frame that drives openms_insight cross-linking.
from core.dataframes.psm_tables import (
    PSM_TABLE_SCHEMA,
    build_psm_table,
)
# The LinePlot/SequenceView/Heatmap LazyFrames keyed off the PSM selected in the Table.
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
