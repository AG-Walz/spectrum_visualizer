"""LazyFrames feeding the openms_insight LinePlot, SequenceView, and Heatmap.
Each builder turns the annotated-spectra tuples produced upstream into a single
typed Polars LazyFrame with one row per peak (peaks, heatmap) or per PSM (sequence). 
"""
import logging
import numpy as np
import pandas as pd
import polars as pl
from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException

from core.spectra.loaded_spectrum import LoadedSpectrum

logger = logging.getLogger(__name__)


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
    Takes (psm_id, annotated_frame) pairs and returns a LazyFrame with one
    row per peak; an empty (but typed) frame is returned when no PSM contributes
    any peaks. The pandas mz/intensity/is_annotated columns are recast
    to the float32/int32 dtypes the openms_insight LinePlot expects.
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
    psm_spectra: list[tuple[int, LoadedSpectrum]],
) -> pl.LazyFrame:
    """Build the per-PSM sequence LazyFrame feeding the openms_insight SequenceView.
    Takes (psm_id, LoadedSpectrum) pairs and returns a LazyFrame with one
    row per PSM matching SEQUENCE_SCHEMA; the peptidoform is converted to
    OpenMS parenthesised-mod notation. An empty (but typed) frame is returned when
    no PSMs are given.
    """
    # No PSMs — return a typed-empty LazyFrame.
    if not psm_spectra:
        return pl.DataFrame(schema=SEQUENCE_SCHEMA).lazy()

    # One row per PSM: sequence in SequenceView notation plus the linking psm_id.
    rows = [
        {
            "sequence": _proforma_to_openms_str(spectrum.peptidoform),
            "precursor_charge": int(spectrum.precursor_charge),
            "psm_id": psm_id,
        }
        for psm_id, spectrum in psm_spectra
    ]

    sequence_frame = pl.DataFrame(rows, schema=SEQUENCE_SCHEMA)
    logger.info(f"Sequence DataFrame: {len(rows)} entries")
    return sequence_frame.lazy()


def build_heatmap_dataframe(spectra: list[LoadedSpectrum]) -> pl.LazyFrame:
    """Stack peaks into an RT × m/z LazyFrame feeding the openms_insight Heatmap.
    Each spectrum contributes one row per peak — its scalar retention time is
    broadcast across every m/z value. Returns a LazyFrame matching
    HEATMAP_SCHEMA, or a typed-empty frame when no spectrum has peaks.
    """
    # Per-spectrum chunks: broadcast the scalar RT across the m/z axis.
    chunks: list[pl.DataFrame] = []
    for spectrum in spectra:
        # Skip empty spectra — keeps the concat below clean.
        if len(spectrum.mz) == 0:
            continue
        chunk = pl.DataFrame({
            "x": np.full(len(spectrum.mz), spectrum.rt, dtype="float32"),
            "y": spectrum.mz.astype("float32"),
            "intensity": spectrum.intensity.astype("float32"),
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
    annotated_queries is the list of (psm_id, LoadedSpectrum, annotated_frame)
    tuples returned by run_streaming_annotation; each builder receives the
    projection it needs. Returns {"peaks", "sequence", "heatmap"} LazyFrames,
    plus a "psm" entry when an already-built psm_frame is supplied.
    """
    # Project the tuples down to what each builder consumes.
    frames = {
        "peaks":    build_peaks_dataframe([(psm_id, df) for psm_id, _, df in annotated_queries]),
        "sequence": build_sequence_dataframe([(psm_id, spec) for psm_id, spec, _ in annotated_queries]),
        "heatmap":  build_heatmap_dataframe([spec for _, spec, _ in annotated_queries]),
    }
    # Pass through a precomputed PSM table when the caller already has one.
    if psm_frame is not None:
        frames["psm"] = psm_frame
    return frames
