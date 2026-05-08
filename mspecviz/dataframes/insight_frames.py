"""Peaks / Sequence / Heatmap LazyFrames for the openms_insight LinePlot, SequenceView, Heatmap."""
import logging

import numpy as np
import pandas as pd
import polars as pl
from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException

from mspecviz.spectra.loaded_spectrum import LoadedSpectrum

logger = logging.getLogger(__name__)


PEAKS_SCHEMA: dict[str, pl.DataType] = {
    "x": pl.Float32,
    "y": pl.Float32,
    "highlight": pl.Int32,
    "annotation": pl.String,
    "ion_type": pl.String,
    "psm_id": pl.Int32,
}

SEQUENCE_SCHEMA: dict[str, pl.DataType] = {
    "sequence": pl.String,
    "precursor_charge": pl.Int32,
    "psm_id": pl.Int32,
}

HEATMAP_SCHEMA: dict[str, pl.DataType] = {
    "x": pl.Float32,
    "y": pl.Float32,
    "intensity": pl.Float32,
}


def _proforma_to_openms_str(peptidoform_str: str) -> str:
    """Convert ProForma (``PEPTM[Oxidation]IDEK``) to OpenMS parens (``PEPTM(Oxidation)IDEK``).

    Feeds SequenceView's ``parse_openms_sequence()``; returns the bare sequence
    on parse failure.
    """
    try:
        peptidoform_obj = Peptidoform(peptidoform_str)
    except PeptidoformException:
        return peptidoform_str.split("/")[0]

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
    """Build a peaks LazyFrame (one row per peak) matching :data:`PEAKS_SCHEMA`.

    ``annotated_spectra`` is ``[(psm_id, frame_from_prepare_spectrum_for_plot)]``.
    ``psm_id`` is the foreign key used as LinePlot/SequenceView ``filters``.
    """
    chunks: list[pl.DataFrame] = []
    for psm_id, annotated_frame in annotated_spectra:
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

    if not chunks:
        return pl.DataFrame(schema=PEAKS_SCHEMA).lazy()

    combined = pl.concat(chunks)
    logger.info(f"Peaks DataFrame: {len(combined)} peaks from {len(chunks)} spectra")
    return combined.lazy()


def build_sequence_dataframe(
    psm_spectra: list[tuple[int, LoadedSpectrum]],
) -> pl.LazyFrame:
    """Build a sequence LazyFrame (one row per PSM) for openms_insight SequenceView.

    ``sequence`` is OpenMS parens notation (e.g. ``SHC(Carbamidomethyl)IAEVEK``)
    so SequenceView's ``parse_openms_sequence()`` renders modifications.
    """
    if not psm_spectra:
        return pl.DataFrame(schema=SEQUENCE_SCHEMA).lazy()

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
    """Build an RT × m/z heatmap LazyFrame (one row per peak) for the openms_insight Heatmap.

    ``spectra`` is not filtered — pass only queries (or both) upstream. Heatmap
    downsamples via multi-resolution streaming internally.
    """
    chunks: list[pl.DataFrame] = []
    for spectrum in spectra:
        if len(spectrum.mz) == 0:
            continue
        chunk = pl.DataFrame({
            "x": np.full(len(spectrum.mz), spectrum.rt, dtype="float32"),
            "y": spectrum.mz.astype("float32"),
            "intensity": spectrum.intensity.astype("float32"),
        })
        chunks.append(chunk)

    if not chunks:
        return pl.DataFrame(schema=HEATMAP_SCHEMA).lazy()

    combined = pl.concat(chunks)
    logger.info(f"Heatmap DataFrame: {len(combined)} points from {len(chunks)} spectra")
    return combined.lazy()
