"""
End-to-end pipeline test for MSInsight.

Exercises peptide search, spectrum loading, peak annotation, and parquet
building using real data files.  This is a plain Python script (no Streamlit).

Usage:
    cd /mnt/volume/nacho/MSInsight
    /mnt/volume/bin/miniconda/envs/msinsight/bin/python test_pipeline.py
"""

import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import polars as pl

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("test_pipeline")


# ---------------------------------------------------------------------------
# Paths to test data
# ---------------------------------------------------------------------------
IDXML_PATH = "/mnt/volume/nacho/data/ProCa_runs/small_dataset/ProT001_1_pout_filtered.idXML"
MZML_DIR = "/mnt/volume/nacho/data/ProCa_runs/small_dataset"
MZML_GLOB = "231204_DF_ProCa_ProT001_*.mzML"
PEPTIDE_LIST_PATH = "/mnt/volume/nacho/data/peptide_lists/peptides.txt"

OUTPUT_DIR = "/mnt/volume/nacho/MSInsight/test_output"

# Cap the number of spectra we process to keep the test fast
MAX_SPECTRA = 50


def main() -> None:
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Validate test data exists
    # ------------------------------------------------------------------
    logger.info("Checking test data paths ...")
    for p in [IDXML_PATH, PEPTIDE_LIST_PATH]:
        if not os.path.isfile(p):
            logger.error(f"File not found: {p}")
            sys.exit(1)

    mzml_paths: List[Path] = sorted(Path(MZML_DIR).glob(MZML_GLOB))
    if not mzml_paths:
        logger.error(f"No mzML files matching {MZML_GLOB} in {MZML_DIR}")
        sys.exit(1)
    logger.info(f"Found {len(mzml_paths)} mzML file(s): {[p.name for p in mzml_paths]}")

    # ------------------------------------------------------------------
    # 2. Load peptide list
    # ------------------------------------------------------------------
    with open(PEPTIDE_LIST_PATH) as fh:
        target_peptides = [line.strip() for line in fh if line.strip()]
    logger.info(f"Loaded {len(target_peptides)} target peptides from {PEPTIDE_LIST_PATH}")

    # ------------------------------------------------------------------
    # 3. Run peptide search
    # ------------------------------------------------------------------
    from src.peptide_search import search_peptides_in_idxmls

    logger.info("Running peptide search ...")
    variant_groups = search_peptides_in_idxmls(
        idxml_paths=[IDXML_PATH],
        target_peptides=target_peptides,
        mzml_paths=[str(p) for p in mzml_paths],
        workers=1,
    )
    logger.info(f"Search returned {len(variant_groups)} variant group(s)")
    total_psms = sum(len(vg.matches) for vg in variant_groups.values())
    logger.info(f"Total PSMs across all variants: {total_psms}")
    for vk, vg in list(variant_groups.items())[:5]:
        logger.info(f"  {vk}: {vg.display_label}  ({len(vg.matches)} PSMs)")

    if total_psms == 0:
        logger.warning("No PSMs found -- nothing to build. Exiting.")
        sys.exit(0)

    # ------------------------------------------------------------------
    # 4. Build PSM table
    # ------------------------------------------------------------------
    from src.parquet_builder import (
        build_psm_table,
        build_peaks_dataframe,
        build_sequence_dataframe,
        build_heatmap_dataframe,
        save_to_parquet,
        validate_dataframe,
        PSM_TABLE_SCHEMA,
        PEAKS_SCHEMA,
        SEQUENCE_SCHEMA,
        HEATMAP_SCHEMA,
    )

    psm_lf, psm_lookup = build_psm_table(variant_groups)
    psm_collected = psm_lf.collect()
    logger.info(f"PSM table shape: {psm_collected.shape}")
    print("\n--- PSM table (first 10 rows) ---")
    print(psm_collected.head(10))
    print()

    # ------------------------------------------------------------------
    # 5. Stream spectra and annotate
    # ------------------------------------------------------------------
    from src.spectrum_loader import load_spectra_streaming, spectrum_lookup_key
    from src.peak_annotation import prepare_spectrum_for_plot

    annotated_list: List[Tuple[int, object, pd.DataFrame]] = []
    spectra_count = 0

    def on_query_ready(result):
        nonlocal spectra_count
        if spectra_count >= MAX_SPECTRA:
            return

        key = spectrum_lookup_key(result.query_spectrum)
        psm_id = psm_lookup.get(key)
        if psm_id is None:
            return

        try:
            annotated_df = prepare_spectrum_for_plot(
                result.query_spectrum.mz,
                result.query_spectrum.intensity,
                result.query_spectrum.peptidoform,
                result.query_spectrum.precursor_charge,
                ion_types=["b", "y"],
                mz_tolerance=0.02,
                normalize=True,
            )
            annotated_list.append((psm_id, result.query_spectrum, annotated_df))
            spectra_count += 1
        except Exception as e:
            logger.warning(f"Annotation failed for psm_id={psm_id}: {e}")

    def on_file_start(filename, is_ref):
        logger.info(f"  Loading file: {filename} (reference={is_ref})")

    def on_file_complete(filename, loaded, not_found):
        logger.info(f"  Completed: {filename} -- loaded={loaded}, not_found={not_found}")

    def on_error(msg):
        logger.error(f"  Loader error: {msg}")

    logger.info("Streaming spectra ...")
    load_spectra_streaming(
        variant_groups=variant_groups,
        mzml_paths=mzml_paths,
        reference_substring="",
        on_query_ready=on_query_ready,
        on_file_start=on_file_start,
        on_file_complete=on_file_complete,
        error_callback=on_error,
        workers=1,
    )
    logger.info(f"Annotated {len(annotated_list)} spectra (cap={MAX_SPECTRA})")

    # ------------------------------------------------------------------
    # 6. Build all four DataFrames
    # ------------------------------------------------------------------
    peaks_lf = build_peaks_dataframe([(pid, df) for pid, _, df in annotated_list])
    seq_lf = build_sequence_dataframe([(pid, s) for pid, s, _ in annotated_list])
    heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated_list])

    # Collect for inspection
    peaks_collected = peaks_lf.collect()
    seq_collected = seq_lf.collect()
    heatmap_collected = heatmap_lf.collect()

    print("\n--- Peaks DataFrame ---")
    print(f"Shape: {peaks_collected.shape}")
    print(f"Schema: {peaks_collected.schema}")
    print(peaks_collected.head(5))

    print("\n--- Sequence DataFrame ---")
    print(f"Shape: {seq_collected.shape}")
    print(f"Schema: {seq_collected.schema}")
    print(seq_collected.head(5))

    print("\n--- Heatmap DataFrame ---")
    print(f"Shape: {heatmap_collected.shape}")
    print(f"Schema: {heatmap_collected.schema}")
    print(heatmap_collected.head(5))

    # ------------------------------------------------------------------
    # 7. Validate schemas
    # ------------------------------------------------------------------
    logger.info("Validating schemas ...")
    all_valid = True
    all_valid &= validate_dataframe(psm_lf, PSM_TABLE_SCHEMA, "PSM table")
    all_valid &= validate_dataframe(peaks_lf, PEAKS_SCHEMA, "Peaks")
    all_valid &= validate_dataframe(seq_lf, SEQUENCE_SCHEMA, "Sequences")
    all_valid &= validate_dataframe(heatmap_lf, HEATMAP_SCHEMA, "Heatmap")
    if all_valid:
        logger.info("All schemas valid.")
    else:
        logger.warning("Some schema mismatches detected (see warnings above).")

    # ------------------------------------------------------------------
    # 8. Save to Parquet
    # ------------------------------------------------------------------
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    save_to_parquet(psm_lf, os.path.join(OUTPUT_DIR, "psm_table.parquet"))
    save_to_parquet(peaks_lf, os.path.join(OUTPUT_DIR, "peaks.parquet"))
    save_to_parquet(seq_lf, os.path.join(OUTPUT_DIR, "sequences.parquet"))
    save_to_parquet(heatmap_lf, os.path.join(OUTPUT_DIR, "heatmap.parquet"))

    # ------------------------------------------------------------------
    # 9. Final summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    print(f"  PSM table  : {psm_collected.shape[0]:>8} rows")
    print(f"  Peaks      : {peaks_collected.shape[0]:>8} rows")
    print(f"  Sequences  : {seq_collected.shape[0]:>8} rows")
    print(f"  Heatmap    : {heatmap_collected.shape[0]:>8} rows")
    print(f"  Elapsed    : {elapsed:>8.1f} s")
    print(f"  Output dir : {OUTPUT_DIR}")
    print("=" * 60)
    logger.info("Pipeline test completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logger.exception("Pipeline test failed with an unhandled exception")
        sys.exit(1)
