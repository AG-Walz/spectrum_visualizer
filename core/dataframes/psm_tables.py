"""Per-PSM LazyFrame + int-keyed lookup for the openms_insight Table and export pipeline (built by ``pipeline.run_streaming_annotation``)."""
import logging

import polars as pl
from psm_utils import PSMList

logger = logging.getLogger(__name__)


PSM_TABLE_SCHEMA: dict[str, pl.DataType] = {
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


def build_psm_table(
    variant_groups: dict[str, PSMList],
) -> tuple[pl.LazyFrame, dict[str, int]]:
    """Build the per-PSM LazyFrame matching :data:`PSM_TABLE_SCHEMA` plus the ``"spectrum_id|peptidoform" -> psm_id`` lookup used to link peaks back to their PSM."""
    rows: list[dict] = []
    psm_lookup: dict[str, int] = {}
    psm_id = 0

    # Walk variants and their PSMs: emit one row per PSM, carrying variant-level fields for grouping
    for variant_key, psm_list in variant_groups.items():
        if not psm_list:
            continue
        base_sequence = psm_list[0].peptidoform.sequence.upper()
        modification = variant_key if variant_key != base_sequence else "(unmodified)"

        for psm in psm_list:
            peptidoform_str = psm.peptidoform.modified_sequence
            # Key must match spectrum_lookup_key() in spectrum_loader
            psm_lookup[f"{psm.spectrum_id}|{peptidoform_str}"] = psm_id
            rows.append({
                "id": psm_id,
                "variant_key": variant_key,
                "base_sequence": base_sequence,
                "modification": modification,
                "display_label": variant_key,
                "peptidoform": peptidoform_str,
                "charge": psm.get_precursor_charge() or 0,
                "score": float(psm.score) if psm.score is not None else 0.0,
                "source_file": psm.run or "",
                "spectrum_id": psm.spectrum_id,
                "idxml_file": psm.metadata.get("idxml_file", ""),
            })
            psm_id += 1

    # Materialise schema-enforced LazyFrame and return alongside the lookup
    df = pl.DataFrame(rows, schema=PSM_TABLE_SCHEMA).lazy()
    logger.info(f"PSM table: {len(rows)} rows across {len(variant_groups)} variants")
    return df, psm_lookup
