"""Peaks / Sequence / Heatmap LazyFrames for openms_insight.

Builds the three visualisation-oriented DataFrames consumed by the LinePlot,
SequenceView, and Heatmap components. Schemas match openms_insight's
``parse_openms_sequence`` / LinePlot expectations exactly.
"""
import logging
import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import polars as pl
from psm_utils import Peptidoform

from mspecviz.dataframes.psm_tables import _enforce_schema
from mspecviz.spectra.loader_models import LoadedSpectrum

logger = logging.getLogger(__name__)


PEAKS_SCHEMA: Dict[str, pl.DataType] = {
    "x": pl.Float32,
    "y": pl.Float32,
    "highlight": pl.Int32,
    "annotation": pl.String,
    "ion_type": pl.String,
    "psm_id": pl.Int32,
}

SEQUENCE_SCHEMA: Dict[str, pl.DataType] = {
    "sequence": pl.String,
    "precursor_charge": pl.Int32,
    "psm_id": pl.Int32,
}

HEATMAP_SCHEMA: Dict[str, pl.DataType] = {
    "x": pl.Float32,
    "y": pl.Float32,
    "intensity": pl.Float32,
}


def _proforma_to_openms_str(peptidoform_str: str) -> str:
    """
    Convert a ProForma peptidoform string to OpenMS parenthesis notation.

    Examples:
        "PEPTIDE"                   -> "PEPTIDE"
        "PEPTM[Oxidation]IDEK"      -> "PEPTM(Oxidation)IDEK"
        "[Acetyl]-PEPTIDE"          -> ".(Acetyl)PEPTIDE"

    Falls back to the base sequence on any parsing error.
    Used to build the ``sequence`` column for openms_insight SequenceView,
    which uses its own ``parse_openms_sequence()`` to handle this format.
    """
    try:
        clean_peptidoform_str = peptidoform_str.split("/")[0]
        peptidoform_obj = Peptidoform(clean_peptidoform_str)
        parts = []
        if peptidoform_obj.properties.get("n_term"):
            for m in peptidoform_obj.properties["n_term"]:
                parts.append(f".({m.name})")
        for aa, residue_mods in peptidoform_obj.parsed_sequence:
            parts.append(aa)
            if residue_mods:
                for m in residue_mods:
                    parts.append(f"({m.name})")
        if peptidoform_obj.properties.get("c_term"):
            for m in peptidoform_obj.properties["c_term"]:
                parts.append(f".({m.name})")
        return "".join(parts)
    except Exception:
        # Safe fallback: return unmodified sequence
        try:
            return Peptidoform(peptidoform_str.split("/")[0]).sequence
        except Exception:
            return peptidoform_str.split("/")[0]


def build_peaks_dataframe(
    annotated_spectra: List[Tuple[int, pd.DataFrame]],
) -> pl.LazyFrame:
    """
    Build a peaks LazyFrame for openms_insight LinePlot and SequenceView.

    Args:
        annotated_spectra: List of (psm_id, annotated_df) where annotated_df is
            the output of mspecviz.annotation.prepare_spectrum_for_plot().
            Expected columns: mz, intensity, annotation, is_annotated, ion_type.

    Returns:
        Polars LazyFrame matching ``PEAKS_SCHEMA`` — one row per peak. Columns:
        x (m/z), y (intensity), highlight (1 if annotated), annotation
        (e.g. ``"y7"``), ion_type (e.g. ``"b"``, ``"precursor"``), and
        psm_id (foreign key into the PSM table — use as the ``filters`` key
        in LinePlot / SequenceView so row-click updates the spectrum view).
    """
    if not annotated_spectra:
        return pl.DataFrame(schema=PEAKS_SCHEMA).lazy()

    chunks: List[pl.DataFrame] = []
    for psm_id, annotated_df in annotated_spectra:
        if annotated_df is None or len(annotated_df) == 0:
            continue
        n = len(annotated_df)

        is_annotated = (
            annotated_df["is_annotated"].values.astype("int32")
            if "is_annotated" in annotated_df.columns
            else np.zeros(n, dtype="int32")
        )
        annotations = (
            annotated_df["annotation"].tolist()
            if "annotation" in annotated_df.columns
            else ["unannotated"] * n
        )
        ion_types = (
            annotated_df["ion_type"].tolist()
            if "ion_type" in annotated_df.columns
            else ["unannotated"] * n
        )

        chunk = pl.DataFrame({
            "x": annotated_df["mz"].values.astype("float32"),
            "y": annotated_df["intensity"].values.astype("float32"),
            "highlight": is_annotated,
            "annotation": annotations,
            "ion_type": ion_types,
            "psm_id": np.full(n, psm_id, dtype="int32"),
        })
        chunks.append(chunk)

    if not chunks:
        return pl.DataFrame(schema=PEAKS_SCHEMA).lazy()

    combined = _enforce_schema(pl.concat(chunks), PEAKS_SCHEMA)
    total_peaks = len(combined)
    logger.info(f"Peaks DataFrame: {total_peaks} peaks from {len(chunks)} spectra")
    return combined.lazy()


def build_sequence_dataframe(
    psm_spectra: List[Tuple[int, "LoadedSpectrum"]],
) -> pl.LazyFrame:
    """
    Build a sequence LazyFrame for the openms_insight SequenceView component.

    One row per PSM. The ``sequence`` column uses OpenMS parenthesis notation
    (e.g. ``SHC(Carbamidomethyl)IAEVEK``) so that openms_insight's
    ``parse_openms_sequence()`` can correctly render modification markers.

    Args:
        psm_spectra: List of (psm_id, LoadedSpectrum) tuples.

    Returns:
        Polars LazyFrame matching ``SEQUENCE_SCHEMA``: sequence (OpenMS
        notation), precursor_charge, psm_id (foreign key). Pass alongside
        ``peaks_lf`` to ``SequenceView`` with ``filters={"selected_psm":
        "psm_id"}`` for linked row-click spectrum rendering.
    """
    if not psm_spectra:
        return pl.DataFrame(schema=SEQUENCE_SCHEMA).lazy()

    rows = [
        {
            "sequence": _proforma_to_openms_str(spec.peptidoform),
            "precursor_charge": int(spec.precursor_charge),
            "psm_id": psm_id,
        }
        for psm_id, spec in psm_spectra
    ]

    df = _enforce_schema(pl.DataFrame(rows), SEQUENCE_SCHEMA)
    logger.info(f"Sequence DataFrame: {len(rows)} entries")
    return df.lazy()


def build_heatmap_dataframe(
    spectra: List["LoadedSpectrum"],
    min_rt: float = 0.0,
) -> pl.LazyFrame:
    """
    Build an RT x m/z heatmap LazyFrame for the openms_insight Heatmap component.

    Each peak from every spectrum contributes one row. The Heatmap component
    handles downsampling internally (multi-resolution streaming).

    Args:
        spectra: List of LoadedSpectrum objects. Reference spectra (is_reference=True)
            are included alongside query spectra -- filter them upstream if needed.
        min_rt: Exclude spectra with RT < min_rt (default 0.0 keeps everything).
            Set to e.g. 60.0 to exclude spectra without a valid RT.

    Returns:
        Polars LazyFrame matching ``HEATMAP_SCHEMA``: x (RT in seconds),
        y (m/z), intensity. One row per peak per spectrum — the Heatmap
        component downsamples internally via multi-resolution streaming.
    """
    chunks: List[pl.DataFrame] = []
    for spec in spectra:
        if spec.rt < min_rt or len(spec.mz) == 0:
            continue
        n = len(spec.mz)
        chunk = pl.DataFrame({
            "x": np.full(n, spec.rt, dtype="float32"),
            "y": spec.mz.astype("float32"),
            "intensity": spec.intensity.astype("float32"),
        })
        chunks.append(chunk)

    if not chunks:
        return pl.DataFrame(schema=HEATMAP_SCHEMA).lazy()

    combined = _enforce_schema(pl.concat(chunks), HEATMAP_SCHEMA)
    total_points = len(combined)
    logger.info(f"Heatmap DataFrame: {total_points} points from {len(chunks)} spectra")
    return combined.lazy()


def save_to_parquet(lf: pl.LazyFrame, path: str) -> None:
    """
    Collect a LazyFrame and write it to a Parquet file.

    Use this when passing ``data_path=`` to openms_insight components instead
    of ``data=``.  Subprocess preprocessing (memory isolation) only works with
    the ``data_path=`` mode.

    Args:
        lf: Polars LazyFrame to persist.
        path: Output file path (e.g. "cache/psm_table.parquet").
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    lf.collect().write_parquet(path, compression="zstd")
    logger.info(f"Saved parquet: {path}")
