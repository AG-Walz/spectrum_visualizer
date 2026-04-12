"""
Build Polars LazyFrames for openms_insight components.

Produces four DataFrames: PSM table, Peaks (LinePlot/SequenceView),
Sequences (SequenceView), and Heatmap (RT x m/z).
"""

import logging
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import polars as pl
from psm_utils import Peptidoform

from src.peptide_search import VariantGroup
from src.spectrum_loader import LoadedSpectrum

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Schema constants  (Dict[str, pl.DataType])
# ---------------------------------------------------------------------------

PSM_TABLE_SCHEMA: Dict[str, pl.DataType] = {
    "id": pl.Int32,
    "variant_key": pl.String,
    "base_sequence": pl.String,
    "modification": pl.String,
    "display_label": pl.String,
    "peptidoform": pl.String,
    "charge": pl.Int32,
    "score": pl.Float32,
    "source_file": pl.String,
    "spectrum_id": pl.String,
    "idxml_file": pl.String,
}

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


# ---------------------------------------------------------------------------
# Schema validation helper
# ---------------------------------------------------------------------------

def validate_dataframe(lf: pl.LazyFrame, schema: dict, name: str) -> bool:
    """Check that a LazyFrame matches the expected schema, logging mismatches."""
    actual_schema = lf.collect_schema()
    valid = True

    for col_name, expected_dtype in schema.items():
        if col_name not in actual_schema:
            logger.warning(f"{name}: missing column '{col_name}' (expected {expected_dtype})")
            valid = False
        elif actual_schema[col_name] != expected_dtype:
            logger.warning(
                f"{name}: column '{col_name}' has type {actual_schema[col_name]}, "
                f"expected {expected_dtype}"
            )
            valid = False

    extra_cols = set(actual_schema.keys()) - set(schema.keys())
    if extra_cols:
        logger.warning(f"{name}: unexpected extra columns: {extra_cols}")
        valid = False

    return valid


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _enforce_schema(df: pl.DataFrame, schema: Dict[str, pl.DataType]) -> pl.DataFrame:
    """Cast all columns to the types declared in *schema*."""
    return df.with_columns([
        pl.col(name).cast(dtype) for name, dtype in schema.items()
        if name in df.columns
    ])


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


# ---------------------------------------------------------------------------
# 1. PSM table  (-> openms_insight Table)
# ---------------------------------------------------------------------------

def build_psm_table(
    variant_groups: Dict[str, VariantGroup],
) -> Tuple[pl.LazyFrame, Dict[str, int]]:
    """
    Build a PSM-level LazyFrame for the openms_insight Table component.

    Assigns a stable integer ``id`` to each PeptideMatch.  The returned
    ``psm_lookup`` dict maps a spectrum key to that ``id`` so that later
    steps (build_peaks_dataframe, build_sequence_dataframe) can link peaks
    back to the correct table row.

    The spectrum key is produced by ``spectrum_lookup_key()`` from
    spectrum_loader:  ``f"{spectrum_id}|{peptidoform}"``.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup from peptide_search.

    Returns:
        Tuple of:
        - LazyFrame with columns:
            id (Int32), variant_key, base_sequence, modification,
            display_label, peptidoform, charge (Int32), score (Float32),
            source_file, spectrum_id, idxml_file
        - psm_lookup dict: spectrum_key -> psm_id (int)
    """
    rows: List[dict] = []
    psm_lookup: Dict[str, int] = {}
    psm_id = 0

    for variant_key, variant_group in variant_groups.items():
        for match in variant_group.matches:
            # Key must match spectrum_lookup_key() in spectrum_loader
            key = f"{match.spectrum_id}|{match.peptidoform}"
            psm_lookup[key] = psm_id
            rows.append({
                "id": psm_id,
                "variant_key": variant_key,
                "base_sequence": variant_group.base_sequence,
                "modification": variant_group.modification_key or "(unmodified)",
                "display_label": variant_group.display_label,
                "peptidoform": match.peptidoform,
                "charge": match.charge,
                "score": float(match.score),
                "source_file": match.mzml_filename or "",
                "spectrum_id": match.spectrum_id,
                "idxml_file": match.idxml_file,
            })
            psm_id += 1

    if not rows:
        return pl.DataFrame(schema=PSM_TABLE_SCHEMA).lazy(), psm_lookup

    df = _enforce_schema(pl.DataFrame(rows), PSM_TABLE_SCHEMA)
    logger.info(f"PSM table: {len(rows)} rows across {len(variant_groups)} variants")
    return df.lazy(), psm_lookup


# ---------------------------------------------------------------------------
# 2. Peaks DataFrame  (-> openms_insight LinePlot / SequenceView peaks_data)
# ---------------------------------------------------------------------------

def build_peaks_dataframe(
    annotated_spectra: List[Tuple[int, pd.DataFrame]],
) -> pl.LazyFrame:
    """
    Build a peaks LazyFrame for openms_insight LinePlot and SequenceView.

    Args:
        annotated_spectra: List of (psm_id, annotated_df) where annotated_df is
            the output of src.peak_annotation.prepare_spectrum_for_plot().
            Expected columns: mz, intensity, annotation, is_annotated, ion_type.

    Returns:
        Polars LazyFrame with columns:
        - x (Float32): m/z value
        - y (Float32): Intensity (0-1 normalized or raw, as prepared upstream)
        - highlight (Int32): 1 if peak is annotated, 0 if unannotated
        - annotation (String): Ion label e.g. "y7", "b3++" or "unannotated"
        - ion_type (String): Ion series e.g. "b", "y", "a", "precursor", "unannotated"
        - psm_id (Int32): Foreign key linking to the PSM table ``id`` column.
            Use this as the ``filters`` key in LinePlot / SequenceView so that
            clicking a row in the Table updates the spectrum view.

    openms_insight LinePlot usage:
        LinePlot(
            cache_id="spectrum",
            data=peaks_lf,
            x_column="x",
            y_column="y",
            highlight_column="highlight",
            annotation_column="annotation",
            filters={"selected_psm": "psm_id"},
            cache_path=CACHE_PATH,
        )
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


# ---------------------------------------------------------------------------
# 3. Sequence DataFrame  (-> openms_insight SequenceView sequence_data)
# ---------------------------------------------------------------------------

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
        Polars LazyFrame with columns:
        - sequence (String): Peptide in OpenMS notation
        - precursor_charge (Int32): Precursor charge state
        - psm_id (Int32): Foreign key for cross-component linking

    openms_insight SequenceView usage (paired with LinePlot):
        sv = SequenceView(
            cache_id="seq_view",
            sequence_data=seq_lf,
            peaks_data=peaks_lf,         # the peaks LazyFrame from build_peaks_dataframe
            filters={"selected_psm": "psm_id"},
            cache_path=CACHE_PATH,
        )
        plot = LinePlot.from_sequence_view(sv, cache_id="spectrum", cache_path=CACHE_PATH)
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


# ---------------------------------------------------------------------------
# 4. Heatmap DataFrame  (-> openms_insight Heatmap)
# ---------------------------------------------------------------------------

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
        Polars LazyFrame with columns:
        - x (Float32): Retention time in seconds
        - y (Float32): m/z value
        - intensity (Float32): Peak intensity

    openms_insight Heatmap usage:
        heatmap = Heatmap(
            cache_id="rt_mz_map",
            data=heatmap_lf,
            x_column="x",
            y_column="y",
            intensity_column="intensity",
            x_label="Retention Time (s)",
            y_label="m/z",
            log_scale=True,
            cache_path=CACHE_PATH,
        )
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


# ---------------------------------------------------------------------------
# Convenience: save DataFrames to Parquet (for openms_insight data_path= mode)
# ---------------------------------------------------------------------------

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
