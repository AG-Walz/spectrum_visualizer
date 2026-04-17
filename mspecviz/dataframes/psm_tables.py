"""PSM-table LazyFrame construction + schema helpers.

Produces the Table-level DataFrame consumed by the openms_insight Table
component, plus the ``psm_lookup`` dict that downstream peaks/sequence builders
use to link peaks to their originating PSM.
"""
import logging
from typing import Dict, List, Tuple

import polars as pl

from mspecviz.search.models import VariantGroup

logger = logging.getLogger(__name__)


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


def validate_dataframe(frame: pl.LazyFrame, schema: dict, name: str) -> bool:
    """Check that a LazyFrame matches the expected schema, logging mismatches."""
    actual_schema = frame.collect_schema()
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


def _enforce_schema(df: pl.DataFrame, schema: Dict[str, pl.DataType]) -> pl.DataFrame:
    """Cast all columns to the types declared in *schema*."""
    return df.with_columns([
        pl.col(name).cast(dtype) for name, dtype in schema.items()
        if name in df.columns
    ])


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
