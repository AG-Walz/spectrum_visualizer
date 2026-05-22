"""Builds the per-PSM table that drives the openms_insight Table selector.
Flattens peptidoform variant_groups into one row per PSMList and produces two outputs: 
a schema-enforced Polars LazyFrame matching PSM_TABLE_SCHEMA 
(selectable PSM list shown in the openms_insight Table) and a "spectrum_id|peptidoform" -> psm_id lookup 
that lets downstream code attach peaks, sequence and heatmap rows back to the integer PSM id that cross-links every component.
(imported by core.dataframes, cli.main and streamlit_app.views.load_annotate.load_annotate_actions.)
"""
import logging
import polars as pl
from psm_utils import PSMList

logger = logging.getLogger(__name__)


# Column layout of the per-PSM LazyFrame — passed verbatim as the Polars schema so
# every build yields identical dtypes regardless of the input rows. id is the
# integer PSM key the openms_insight components cross-link on.
PSM_TABLE_SCHEMA: dict[str, pl.DataType] = {
    "id": pl.Int32,             # Sequential 0-based PSM id — the cross-link key
    "variant_key": pl.String,   # Peptidoform variant this PSM belongs to (group label)
    "base_sequence": pl.String, # Bare amino-acid sequence, uppercase, no modifications
    "modification": pl.String,  # Modification descriptor, or "(unmodified)"
    "display_label": pl.String, # Human-facing label shown in the Table (mirrors variant_key)
    "peptidoform": pl.String,   # Full modified sequence (ProForma string)
    "charge": pl.Int32,         # Precursor charge state
    "score": pl.Float32,        # Search-engine score for this PSM
    "source_file": pl.String,   # mzML run the spectrum came from
    "spectrum_id": pl.String,   # Native spectrum identifier within the run
    "idxml_file": pl.String,    # idXML file the PSM was read from
}


def build_psm_table(
    variant_groups: dict[str, PSMList],
) -> tuple[pl.LazyFrame, dict[str, int]]:
    """Flatten variant_groups into the per-PSM table and its id lookup. Each key in
    variant_groups maps each peptidoform variant key to a PSMList; every PSM
    becomes one row carrying both PSM-level fields and the variant-level grouping
    fields shared by its siblings. Returns (lazyframe, psm_lookup) where the
    LazyFrame matches PSM_TABLE_SCHEMA and psm_lookup maps
    "spectrum_id|peptidoform" to the integer id — the same key shape that
    spectrum_lookup_key() produces, so loaded spectra resolve to their PSM row.
    """
    rows: list[dict] = []
    psm_lookup: dict[str, int] = {}
    # PSM ids are assigned sequentially across all variants — global, not per-variant.
    psm_id = 0

    # Walk variants and their PSMs: emit one row per PSM, carrying variant-level fields for grouping
    for variant_key, psm_list in variant_groups.items():
        if not psm_list:
            continue
        # Variant-level fields — derived once and repeated on every PSM row of the group.
        # A variant_key equal to the bare sequence means no modifications were applied.
        base_sequence = psm_list[0].peptidoform.sequence.upper()
        modification = variant_key if variant_key != base_sequence else "(unmodified)"

        for psm in psm_list:
            peptidoform_str = psm.peptidoform.modified_sequence
            # Key must match spectrum_lookup_key() in spectrum_loader
            psm_lookup[f"{psm.spectrum_id}|{peptidoform_str}"] = psm_id
            # Coerce missing/optional PSM attributes to schema-safe defaults so the
            # Polars schema never has to widen a column to Null.
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
