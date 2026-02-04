#!/usr/bin/env python3
"""
MSpecViz Pipeline Runner for Profiling

Run with cProfile:
    python -m cProfile -o program.prof main.py

Visualize with SnakeViz:
    snakeviz program.prof

"""
import time
import json
from pathlib import Path
from datetime import datetime

# ============================================================================
# CONFIGURATION - Edit these values for your test
# ============================================================================

# Dataset paths
IDXML_PATH = Path("/mnt/volume/nacho/data/ProCa_Anna/small_dataset/")
MZML_PATH = Path("/mnt/volume/nacho/data/ProCa_Anna/small_dataset/")

# Target peptides to search
TARGET_PEPTIDES = [
    "KPFDQTTISL",
    "ATELTHEDY",
    "AYIKNGQLF",
    "ATYSVGMQK",
]

# Reference file substring (files containing this are reference spectra)
# Note: matching is case-insensitive
REFERENCE_SUBSTRING = "synpep"

# Analysis settings
ION_TYPES = ["b", "y"]
MZ_TOLERANCE = 0.02  # Da

# Export settings
EXPORT_DIR = Path("/mnt/volume/nacho/MSpecViz/exports/")
EXPORT_FORMATS = ["png", "svg"]

# ============================================================================
# TIMING UTILITIES
# ============================================================================

class Timer:
    """Context manager for timing code sections."""

    def __init__(self, name: str):
        self.name = name
        self.start = None
        self.elapsed = None

    def __enter__(self):
        self.start = time.perf_counter()
        print(f"\n{'='*60}")
        print(f"[START] {self.name}")
        print(f"{'='*60}")
        return self

    def __exit__(self, *args):
        self.elapsed = time.perf_counter() - self.start
        print(f"[DONE] {self.name}: {self.elapsed:.2f}s")
        return False


timings = {}

def record_timing(name: str, elapsed: float):
    """Record timing for summary."""
    timings[name] = elapsed


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Run the complete MSpecViz pipeline."""

    pipeline_start = time.perf_counter()

    print("\n" + "="*70)
    print("MSpecViz Pipeline - Profiling Run")
    print(f"Started: {datetime.now().isoformat()}")
    print("="*70)

    # -------------------------------------------------------------------------
    # Step 0: Import modules (measure import time)
    # -------------------------------------------------------------------------
    with Timer("Import modules") as t:
        from core.peptide_search import search_peptides_in_idxmls, summarize_search_results
        from core.spectrum_loader import load_matched_spectra
        from core.peak_annotation import prepare_spectrum_for_plot
        from core.mirror_plot import create_mirror_plot, figure_to_bytes
        from core.export import export_mirror_plots
    record_timing("imports", t.elapsed)

    # -------------------------------------------------------------------------
    # Step 1: Discover files
    # -------------------------------------------------------------------------
    with Timer("Discover files") as t:
        idxml_files = sorted(IDXML_PATH.glob("*.idXML"))
        mzml_files = sorted(MZML_PATH.glob("*.mzML"))

        print(f"Found {len(idxml_files)} idXML files:")
        for f in idxml_files:
            print(f"  - {f.name}")

        print(f"\nFound {len(mzml_files)} mzML files:")
        for f in mzml_files:
            size_gb = f.stat().st_size / (1024**3)
            print(f"  - {f.name} ({size_gb:.2f} GB)")
    record_timing("file_discovery", t.elapsed)

    # -------------------------------------------------------------------------
    # Step 2: Search peptides in idXML files
    # -------------------------------------------------------------------------
    with Timer("Search peptides in idXML") as t:
        idxml_paths = [str(f) for f in idxml_files]

        variant_groups = search_peptides_in_idxmls(
            idxml_paths=idxml_paths,
            target_peptides=TARGET_PEPTIDES,
        )

        print(f"\nSearch results:")
        print(f"  Target peptides: {len(TARGET_PEPTIDES)}")
        print(f"  Variant groups found: {len(variant_groups)}")

        total_matches = sum(len(vg.matches) for vg in variant_groups.values())
        print(f"  Total PSMs: {total_matches}")

        # Show summary
        if variant_groups:
            summary_df = summarize_search_results(variant_groups)
            print("\nVariant summary:")
            print(summary_df.to_string(index=False))
    record_timing("peptide_search", t.elapsed)

    if not variant_groups:
        print("\nNo peptides found! Exiting.")
        return

    # -------------------------------------------------------------------------
    # Step 3: Load spectra from mzML files
    # -------------------------------------------------------------------------
    with Timer("Load spectra from mzML") as t:
        mzml_paths = list(mzml_files)  # Keep as Path objects

        def progress_callback(progress):
            """Print progress updates."""
            print(f"  Progress: {progress.files_done}/{progress.total_files} files, "
                  f"{progress.spectra_loaded}/{progress.spectra_expected} spectra, "
                  f"{progress.spectra_not_found} not found - {progress.current_file}", end='\r')

        loaded_variants = load_matched_spectra(
            variant_groups=variant_groups,
            mzml_paths=mzml_paths,
            reference_substring=REFERENCE_SUBSTRING,
            progress_callback=progress_callback,
        )

        print()  # New line after progress

        # Summary
        total_query = sum(len(vg.query_spectra) for vg in loaded_variants.values())
        total_ref = sum(len(vg.reference_spectra) for vg in loaded_variants.values())
        print(f"\nLoaded spectra:")
        print(f"  Query spectra: {total_query}")
        print(f"  Reference spectra: {total_ref}")
    record_timing("spectrum_loading", t.elapsed)

    # -------------------------------------------------------------------------
    # Step 4: Prepare spectra for plotting (annotation)
    # -------------------------------------------------------------------------
    with Timer("Annotate spectra") as t:
        annotation_count = 0

        for variant_key, vg in loaded_variants.items():
            # Annotate query spectra
            for spec in vg.query_spectra:
                if spec.peptidoform and spec.mz is not None:
                    try:
                        df = prepare_spectrum_for_plot(
                            mz_array=spec.mz,
                            intensity_array=spec.intensity,
                            peptidoform=spec.peptidoform,
                            charge=spec.precursor_charge,
                            ion_types=ION_TYPES,
                            mz_tolerance=MZ_TOLERANCE,
                        )
                        spec.annotated_df = df
                        annotation_count += 1
                    except Exception as e:
                        print(f"  Warning: Failed to annotate {spec.peptidoform}: {e}")

            # Annotate reference spectra
            for spec in vg.reference_spectra:
                if spec.peptidoform and spec.mz is not None:
                    try:
                        df = prepare_spectrum_for_plot(
                            mz_array=spec.mz,
                            intensity_array=spec.intensity,
                            peptidoform=spec.peptidoform,
                            charge=spec.precursor_charge,
                            ion_types=ION_TYPES,
                            mz_tolerance=MZ_TOLERANCE,
                        )
                        spec.annotated_df = df
                        annotation_count += 1
                    except Exception as e:
                        print(f"  Warning: Failed to annotate {spec.peptidoform}: {e}")

        print(f"Annotated {annotation_count} spectra")
    record_timing("annotation", t.elapsed)

    # -------------------------------------------------------------------------
    # Step 5: Generate mirror plots
    # -------------------------------------------------------------------------
    with Timer("Generate mirror plots") as t:
        plot_count = 0
        figures = []

        for variant_key, vg in loaded_variants.items():
            # Get annotated spectra
            query_annotated = [s for s in vg.query_spectra if hasattr(s, 'annotated_df')]
            ref_annotated = [s for s in vg.reference_spectra if hasattr(s, 'annotated_df')]

            # Generate comparison pairs
            if ref_annotated:
                # Query vs Reference
                for query_spec in query_annotated:
                    for ref_spec in ref_annotated:
                        try:
                            fig = create_mirror_plot(
                                upper_df=query_spec.annotated_df,
                                lower_df=ref_spec.annotated_df,
                                upper_label=f"Query: {Path(query_spec.source_file).stem}",
                                lower_label=f"Ref: {Path(ref_spec.source_file).stem}",
                                title=f"{vg.base_sequence}",
                            )
                            figures.append((variant_key, query_spec, ref_spec, fig))
                            plot_count += 1
                        except Exception as e:
                            print(f"  Warning: Failed to create plot: {e}")
            else:
                # Query vs Query (no reference available)
                for i, spec1 in enumerate(query_annotated):
                    for spec2 in query_annotated[i+1:]:
                        try:
                            fig = create_mirror_plot(
                                upper_df=spec1.annotated_df,
                                lower_df=spec2.annotated_df,
                                upper_label=f"{Path(spec1.source_file).stem}",
                                lower_label=f"{Path(spec2.source_file).stem}",
                                title=f"{vg.base_sequence}",
                            )
                            figures.append((variant_key, spec1, spec2, fig))
                            plot_count += 1
                        except Exception as e:
                            print(f"  Warning: Failed to create plot: {e}")

        print(f"Generated {plot_count} mirror plots")
    record_timing("plot_generation", t.elapsed)

    # -------------------------------------------------------------------------
    # Step 6: Export plots
    # -------------------------------------------------------------------------
    with Timer("Export plots") as t:
        # Ensure export directory exists
        EXPORT_DIR.mkdir(parents=True, exist_ok=True)

        # Export mirror plots (creates timestamped subdirectory automatically)
        exported_paths = export_mirror_plots(
            variant_groups=loaded_variants,
            output_dir=str(EXPORT_DIR),
            formats=EXPORT_FORMATS,
            ion_types=ION_TYPES,
            mz_tolerance=MZ_TOLERANCE,
        )

        total_exported = sum(len(paths) for paths in exported_paths.values())
        print(f"\nExported {total_exported} files:")
        for fmt, paths in exported_paths.items():
            print(f"  {fmt}: {len(paths)} files")

        if exported_paths.get('png'):
            print(f"Export directory: {Path(exported_paths['png'][0]).parent.parent}")
    record_timing("export", t.elapsed)

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    pipeline_elapsed = time.perf_counter() - pipeline_start

    print("\n" + "="*70)
    print("TIMING SUMMARY")
    print("="*70)

    for step, elapsed in timings.items():
        pct = (elapsed / pipeline_elapsed) * 100
        bar = "█" * int(pct / 2)
        print(f"{step:25s}: {elapsed:8.2f}s ({pct:5.1f}%) {bar}")

    print("-"*70)
    print(f"{'TOTAL':25s}: {pipeline_elapsed:8.2f}s")
    print("="*70)

    # Save timing results
    timing_file = EXPORT_DIR / f"timing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(timing_file, 'w') as f:
        json.dump({
            "timings": timings,
            "total": pipeline_elapsed,
            "config": {
                "idxml_path": str(IDXML_PATH),
                "mzml_path": str(MZML_PATH),
                "target_peptides": TARGET_PEPTIDES,
                "reference_substring": REFERENCE_SUBSTRING,
                "ion_types": ION_TYPES,
                "mz_tolerance": MZ_TOLERANCE,
            }
        }, f, indent=2)
    print(f"\nTiming saved to: {timing_file}")

    print("\nDone! Run with cProfile for detailed profiling:")
    print("  python -m cProfile -o program.prof main.py")
    print("  snakeviz program.prof")


if __name__ == "__main__":
    main()
