"""Polars LazyFrames feeding openms_insight Table/LinePlot/SequenceView/Heatmap.
- `build_psm_table` flattens peptidoform_psms into the per-PSM table
  that drives the openms_insight Table selector. It returns a "spectrum_id|peptidoform" 
  psm_id lookup to attach peaks, sequence and heatmap rows back to the integer PSM id.
- `build_peaks_dataframe`, `build_sequence_dataframe`, `build_heatmap_dataframe`
  turn the annotated-spectra tuples produced upstream into the per-peak /
  per-PSM LazyFrames that the LinePlot, SequenceView and Heatmap consume.
- `build_visualization_frames` bundles those three into the dict openms_insight
  expects (and passes a precomputed PSM table through unchanged).

Imported by cli.main and streamlit_app.views.load_annotate.load_annotate_actions.
"""
import logging

import numpy as np
import pandas as pd
import polars as pl
from psm_utils import Peptidoform, PSM, PSMList
from psm_utils.peptidoform import PeptidoformException

from core.spectra.streaming_types import psm_lookup_key

logger = logging.getLogger(__name__)


# ============================================================================
# Schemas
# ============================================================================

# Column layout of the per-PSM LazyFrame — passed verbatim as the Polars schema so
# every build yields identical dtypes regardless of the input rows. id is the
# integer PSM key the openms_insight components cross-link on.
PSM_TABLE_SCHEMA: dict[str, pl.DataType] = {
    "id": pl.Int32,             # Sequential 0-based PSM id — the cross-link key
    "base_sequence": pl.String, # Bare amino-acid sequence, uppercase, no modifications
    "modification": pl.String,  # Modification descriptor, or "(unmodified)"
    "display_label": pl.String, # Human-facing label shown in the Table (mirrors peptidoform)
    "peptidoform": pl.String,   # Full modified sequence (ProForma string) — the grouping key
    "charge": pl.Int32,         # Precursor charge state
    "score": pl.Float32,        # Search-engine score for this PSM
    "source_file": pl.String,   # mzML run the spectrum came from
    "spectrum_id": pl.String,   # Native spectrum identifier within the run
    "idxml_file": pl.String,    # idXML file the PSM was read from
    "source": pl.String,        # Source folder of the idXML — disambiguates same-basename files
}

# LazyFrame schema for the per-peak frame consumed by the openms_insight LinePlot.
PEAKS_SCHEMA: dict[str, pl.DataType] = {
    "x": pl.Float32,          # peak m/z
    "y": pl.Float32,          # peak intensity
    "highlight": pl.Int32,    # 1 when the peak matched a theoretical fragment, else 0
    "annotation": pl.String,  # fragment-ion label (e.g. "b2", "y5")
    "ion_type": pl.String,    # ion series ("b", "y", ...) driving the LinePlot colour
    "psm_id": pl.Int32,       # PSM this peak belongs to — the linking key
}

# LazyFrame schema for the per-PSM frame consumed by the openms_insight SequenceView.
SEQUENCE_SCHEMA: dict[str, pl.DataType] = {
    "sequence": pl.String,           # peptide in OpenMS parenthesised-mod notation
    "precursor_charge": pl.Int32,    # precursor charge state
    "psm_id": pl.Int32,              # PSM this row belongs to — the linking key
}

# LazyFrame schema for the RT × m/z frame consumed by the openms_insight Heatmap.
HEATMAP_SCHEMA: dict[str, pl.DataType] = {
    "x": pl.Float32,          # retention time
    "y": pl.Float32,          # peak m/z
    "intensity": pl.Float32,  # peak intensity
}


# ============================================================================
# PSM selector table
# ============================================================================

def build_psm_table(
    peptidoform_psms: dict[str, PSMList],
) -> tuple[pl.LazyFrame, dict[str, int]]:
    """Flatten peptidoform_psms into the per-PSM table and its id lookup.
    Every PSM becomes one row that carries both its own PSM-level fields and the
    peptidoform-level fields shared by its siblings, matching PSM_TABLE_SCHEMA. Returns
    (lazyframe, psm_lookup), where psm_lookup maps each psm_lookup_key()
    ("spectrum_id|peptidoform") to its integer id, so loaded spectra resolve back to their row.
    """
    rows: list[dict] = []
    psm_lookup: dict[str, int] = {}

    # Walk every PSM under every peptidoform: id is the row index — globally
    # sequential. Each row carries the PSM's own fields plus the peptidoform-level
    # fields shared by its siblings. A peptidoform equal to the bare sequence means
    # unmodified. Optional fields are coerced to schema-safe defaults so the Polars
    # schema never has to widen a column to Null.
    for peptidoform, psm_list in peptidoform_psms.items():
        for psm in psm_list:
            psm_id = len(rows)
            peptidoform_str = psm.peptidoform.modified_sequence
            base_sequence = psm.peptidoform.sequence.upper()
            psm_lookup[psm_lookup_key(psm)] = psm_id
            rows.append({
                "id": psm_id,
                "base_sequence": base_sequence,
                "modification": peptidoform if psm.peptidoform.is_modified else "(unmodified)",
                "display_label": peptidoform,
                "peptidoform": peptidoform_str,
                "charge": psm.get_precursor_charge() or 0,
                "score": float(psm.score) if psm.score is not None else 0.0,
                "source_file": psm.run or "",
                "spectrum_id": psm.spectrum_id,
                "idxml_file": psm.metadata.get("idxml_file", ""),
                "source": psm.metadata.get("idxml_source", ""),
            })

    df = pl.DataFrame(rows, schema=PSM_TABLE_SCHEMA).lazy()
    logger.info(f"PSM table: {len(rows)} rows across {len(peptidoform_psms)} peptidoforms")
    return df, psm_lookup


# ============================================================================
# LinePlot / SequenceView / Heatmap frames
# ============================================================================

def _proforma_to_openms_str(peptidoform_str: str) -> str:
    """Convert a ProForma peptidoform to the parenthesised-mod notation the SequenceView expects.
    (PEPTM[Oxidation]IDEK -> PEPTM(Oxidation)IDEK); When the input cannot be parsed as
    ProForma, the bare sequence is returned unchanged as fallback.
    """
    # Unparseable input — drop any "/charge" suffix and return the raw sequence.
    try:
        peptidoform_obj = Peptidoform(peptidoform_str)
    except PeptidoformException:
        return peptidoform_str.split("/")[0]

    # Emit N-term mods, then per-residue letters with their mods, then C-term mods
    parts = []
    if peptidoform_obj.properties.get("n_term"):
        for mod in peptidoform_obj.properties["n_term"]:
            parts.append(f".({mod.name})")
    for amino_acid, residue_mods in peptidoform_obj.parsed_sequence:
        parts.append(amino_acid)
        if residue_mods:
            for mod in residue_mods:
                parts.append(f"({mod.name})")
    if peptidoform_obj.properties.get("c_term"):
        for mod in peptidoform_obj.properties["c_term"]:
            parts.append(f".({mod.name})")
    return "".join(parts)


def build_peaks_dataframe(
    annotated_spectra: list[tuple[int, pd.DataFrame]],
) -> pl.LazyFrame:
    """Stack annotated peaks across PSMs into one LazyFrame matching PEAKS_SCHEMA.
    Takes (psm_id, annotated_frame) pairs and returns one row per peak — empty
    (but typed) when no PSM contributes peaks. The pandas mz/intensity/is_annotated
    columns are recast to the float32/int32 dtypes the LinePlot expects.
    """
    # Per-PSM chunks: each annotated frame becomes a typed Polars block.
    chunks: list[pl.DataFrame] = []
    for psm_id, annotated_frame in annotated_spectra:
        # Skip PSMs that produced no peaks — keeps the concat below clean.
        if annotated_frame is None or len(annotated_frame) == 0:
            continue
        chunk = pl.DataFrame({
            "x": annotated_frame["mz"].values.astype("float32"),
            "y": annotated_frame["intensity"].values.astype("float32"),
            "highlight": annotated_frame["is_annotated"].values.astype("int32"),
            "annotation": annotated_frame["annotation"].tolist(),
            "ion_type": annotated_frame["ion_type"].tolist(),
            "psm_id": np.full(len(annotated_frame), psm_id, dtype="int32"),
        })
        chunks.append(chunk)

    # No PSM contributed peaks — return a typed-empty LazyFrame so downstream
    # components still see the expected schema.
    if not chunks:
        return pl.DataFrame(schema=PEAKS_SCHEMA).lazy()

    # Concatenate every per-PSM block into the final LazyFrame.
    combined = pl.concat(chunks)
    logger.info(f"Peaks DataFrame: {len(combined)} peaks from {len(chunks)} spectra")
    return combined.lazy()


def build_sequence_dataframe(
    psm_entries: list[tuple[int, PSM]],
) -> pl.LazyFrame:
    """Build the per-PSM sequence LazyFrame feeding the openms_insight SequenceView.
    Takes (psm_id, PSM) pairs and returns one row per PSM matching SEQUENCE_SCHEMA, with the
    peptidoform converted to OpenMS parenthesised-mod notation and the charge taken from the
    peptidoform's identification charge. An empty (but typed) frame is returned when no PSMs are given.
    """
    # No PSMs — return a typed-empty LazyFrame.
    if not psm_entries:
        return pl.DataFrame(schema=SEQUENCE_SCHEMA).lazy()

    # One row per PSM: sequence in SequenceView notation plus the linking psm_id.
    rows = [
        {
            "sequence": _proforma_to_openms_str(psm.peptidoform.modified_sequence),
            # identification charge carried by the peptidoform (= precursor charge for confident PSMs)
            "precursor_charge": int(psm.get_precursor_charge() or 0),
            "psm_id": psm_id,
        }
        for psm_id, psm in psm_entries
    ]

    sequence_frame = pl.DataFrame(rows, schema=SEQUENCE_SCHEMA)
    logger.info(f"Sequence DataFrame: {len(rows)} entries")
    return sequence_frame.lazy()


def build_heatmap_dataframe(psms: list[PSM]) -> pl.LazyFrame:
    """Stack peaks into an RT x m/z LazyFrame feeding the openms_insight Heatmap.
    Each PSM contributes one row per peak (read off psm.spectrum) — its scalar
    retention time is broadcast across every m/z value. Returns a LazyFrame
    matching HEATMAP_SCHEMA, or a typed-empty frame when no PSM has peaks.
    """
    # Per-PSM chunks: broadcast the scalar RT across the m/z axis.
    chunks: list[pl.DataFrame] = []
    for psm in psms:
        mz = psm.spectrum["mz"]
        # Skip empty spectra — keeps the concat below clean.
        if len(mz) == 0:
            continue
        chunk = pl.DataFrame({
            "x": np.full(len(mz), psm.retention_time or 0.0, dtype="float32"),
            "y": mz.astype("float32"),
            "intensity": psm.spectrum["intensity"].astype("float32"),
        })
        chunks.append(chunk)

    # No spectrum had peaks — return a typed-empty LazyFrame.
    if not chunks:
        return pl.DataFrame(schema=HEATMAP_SCHEMA).lazy()

    # Concatenate every per-spectrum block into the final LazyFrame.
    combined = pl.concat(chunks)
    logger.info(f"Heatmap DataFrame: {len(combined)} points from {len(chunks)} spectra")
    return combined.lazy()


def build_visualization_frames(annotated_queries, psm_frame=None) -> dict:
    """Bundle the peaks / sequence / heatmap LazyFrames into the dict openms_insight expects.
    annotated_queries is the (psm_id, query_psm, annotated_frame) list from run_streaming_annotation;
    each builder receives just the projection it needs. Returns {"peaks", "sequence", "heatmap"}
    LazyFrames, plus a "psm" entry when a precomputed psm_frame is supplied.
    """
    # Project the tuples down to what each builder consumes.
    frames = {
        "peaks":    build_peaks_dataframe([(psm_id, df) for psm_id, _, df in annotated_queries]),
        "sequence": build_sequence_dataframe([(psm_id, psm) for psm_id, psm, _ in annotated_queries]),
        "heatmap":  build_heatmap_dataframe([psm for _, psm, _ in annotated_queries]),
    }
    # Pass through a precomputed PSM table when the caller already has one.
    if psm_frame is not None:
        frames["psm"] = psm_frame
    return frames
