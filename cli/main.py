"""MSpecViz CLI entry: run the 3-step pipeline (lookup → load+annotate → export) delegating functionality to core/ packages;
configured from a YAML config invoked via python -m cli <config.yaml>.
"""
import glob as globmod
import logging
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
import click
import pandas as pd
import yaml
from tqdm import tqdm

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Imports from core package
from core.defaults import DEFAULTS  # noqa: E402
from core.dataframes import build_psm_table, build_visualization_frames  # noqa: E402
from core.export import (  # noqa: E402
    ExportCallbacks, ExportSelection, PlotSettings,
    build_run_settings, run_export,
)
from core.pipeline import (  # noqa: E402
    build_annotation_kwargs, build_koina_kwargs, run_streaming_annotation,
)
from core.peptide_lookup import load_idxmls, lookup_peptides  # noqa: E402

logger = logging.getLogger("core.cli")

# CLI-only fields not in DEFAULTS — populated by the user's YAML.
_CLI_EXTRAS: dict = {
    "idxml_files": [], "mzml_files": [],
    "peptide_source": "", "peptides": [],
    "reference_substring": "", "no_plots": False,
}


def _expand_globs(patterns: list) -> list:
    """Resolves shell glob patterns (e.g. "runs/**/ProT*.idXML")  into a sorted, 
    unique (de-duplicated) list of absolute paths. Necessary because YAML lets the user 
    write paths with wildcards — core/ expects concrete file lists.
    """
    result, seen = [], set()
    for pattern in patterns:
        # Sort the YAML-provided patterns for deterministic output 
        for match in sorted(globmod.glob(pattern, recursive=True)):
            # Convert matches to absolute paths 
            absolute = str(Path(match).resolve())
            # De-duplicate across patterns 
            if absolute not in seen:
                seen.add(absolute)
                result.append(absolute)
    return result


def _resolve_peptide_source(source) -> list:
    """Resolve the peptide loading source (a YAML list, or a file path or inline
    string holding one peptide per line) into a list of upper-cased sequences.
    """
    # A. Sequence list case: strip and uppercase each entry
    if isinstance(source, list):
        return [str(peptide).strip().upper() for peptide in source if str(peptide).strip()]
    # B. File path or inline string case: one peptide per line
    source = str(source).strip()
    text = Path(source).read_text() if Path(source).is_file() else source
    return [line.strip().upper() for line in text.splitlines() if line.strip()]


def _load_config(path: str) -> SimpleNamespace:
    """Full config-loading pipeline: Parse the YAML at path, join with the defaults, 
    expand file globs, and resolve the peptide source into a SimpleNamespace.
    """
    # Parse YAML file: error out on non-mapping payload
    config_path = Path(path)
    yaml_config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(yaml_config, dict):
        sys.exit(f"Error: config must be a YAML mapping: {path}")

    # Merge 3 layers: CLI extras (renamed defaults) < DEFAULTS < user YAML (last write wins)
    merged: dict = {**_CLI_EXTRAS, **asdict(DEFAULTS), **yaml_config}
    # Wrap merged dictionary in a SimpleNamespace for convenient attribute access
    config = SimpleNamespace(**merged)

    # Expand input globs and require non-empty idXML and mzML lists
    config.idxml_files = _expand_globs(config.idxml_files)
    config.mzml_files = _expand_globs(config.mzml_files)
    if not config.idxml_files or not config.mzml_files:
        sys.exit("Error: idxml_files and mzml_files must both be set.")
    if not config.peptide_source:
        sys.exit("Error: 'peptide_source' is required.")
    
    # Resolve the peptide source into a list of sequences 
    config.peptides = _resolve_peptide_source(config.peptide_source)
    if not config.peptides:
        sys.exit("Error: peptide list is empty.")
    return config


def _run_lookup(config) -> tuple:
    """Load idXML files once into a PSMList (core.peptide_lookup.load_idxmls), then
    scan that in-memory list for the target peptides (core.peptide_lookup.lookup_peptides),
    producing grouped_psms (PSMList per peptidoform). Returns (grouped_psms, total_psms).
    """
    logger.info("Step 1: Looking up peptides in idXML files")
    # Progress bar
    bar = tqdm(total=len(config.idxml_files), desc="Loading", unit="file")
    errors: list = []
    # Parse every idXML once, then filter/group target peptides from memory (no re-read)
    psm_list = load_idxmls(
        idxml_paths=config.idxml_files,
        progress_callback=lambda current, total, message: bar.update(current - bar.n),
        error_callback=errors.append,
        workers=min(config.workers, len(config.idxml_files)),
    )
    bar.close()
    grouped_psms = lookup_peptides(
        psm_list, target_peptides=config.peptides, mzml_paths=config.mzml_files,
    )
    # Log any lookup errors but keep going
    for error in errors:
        logger.warning(error)
    # Count total PSMs across all grouped PSMs for progress tracking in the next step
    total_psms = sum(len(psm_list) for psm_list in grouped_psms.values())
    logger.info(f"Found {total_psms} PSMs across {len(grouped_psms)} peptidoforms")
    return grouped_psms, total_psms


def _run_load(config, grouped_psms, psm_lookup, total_psms) -> tuple:
    """Calls core.pipeline.run_streaming_annotation() to stream spectra from mzMLs and 
    annotate fragment ions (optionally fetching Koina AI prediction when no reference spectra are available);
    return (annotated_queries, reference_store, comparison_rows).
    """
    logger.info("Step 2: Loading and annotating spectra")
    # Build a simple getter function for accessing config values with defaults
    getter = lambda key, default: getattr(config, key, default)
    # Progress bar
    bar = tqdm(total=total_psms, desc="Loading", unit="spec")
    annotated_queries, reference_store, comparison_rows = run_streaming_annotation(
        grouped_psms=grouped_psms,
        mzml_paths=config.mzml_files,
        reference_substring=config.reference_substring,
        psm_lookup=psm_lookup,
        annotation_kwargs=build_annotation_kwargs(getter),
        koina_kwargs=build_koina_kwargs(getter),
        workers=config.workers,
        on_progress=lambda count: bar.update(count - bar.n),
        on_file_start=lambda path, is_ref: logger.info(
            f"Processing {'reference' if is_ref else 'query'}: {Path(path).name}"),
        on_file_complete=lambda path, loaded, total: logger.info(
            f"Completed {Path(path).name}: {loaded}/{total} spectra"),
        on_error=lambda error: logger.error(f"Streaming error: {error}"),
    )
    bar.close()
    return annotated_queries, reference_store, comparison_rows


def _run_export(config, annotated_queries, reference_store, grouped_psms,
                comparison_rows, psm_frame):
    """Calls core.export.run_export to write parquet/CSV tables and 
    (unless --no-plots is selected) mirror & summary plots.
    """
    logger.info("Step 3: Exporting results")
    # Progress bar 
    plot_bar: dict = {"bar": None}
    def _on_plot_progress(completed: int, total: int) -> None:
        if plot_bar["bar"] is None:
            plot_bar["bar"] = tqdm(total=total, desc="Plots", unit="plot")
        plot_bar["bar"].update(completed - plot_bar["bar"].n)
        if completed >= total:
            plot_bar["bar"].close()

    return run_export(
        export_directory=config.export_directory,
        project_name=config.project_name,
        run_settings=build_run_settings(
            lambda key, default: getattr(config, key, default),
            peptides=config.peptides,
        ),
        frames=build_visualization_frames(annotated_queries, psm_frame),
        annotated_queries=annotated_queries,
        reference_store=reference_store,
        grouped_psms=grouped_psms,
        comparison_frame=pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame(),
        # Export selection: which tables and plots to write based on config
        selection=ExportSelection(
            parquet=True, csv=True,
            mirror_plots=not config.no_plots,
            summary_plots=not config.no_plots,
            html="html" in config.export_formats,
        ),
        # Plot settings: mostly from config (or DEFAULTS)
        plot_settings=PlotSettings(
            export_formats=list(config.export_formats),
            figure_size=tuple(config.figure_size),
            export_resolution=config.export_resolution,
            plot_unpaired_queries=config.plot_unpaired_queries,
            workers=config.workers,
        ),
        # Progress and error reporting callbacks for the export step
        callbacks=ExportCallbacks(
            on_step=logger.info, on_plot_progress=_on_plot_progress,
            on_warning=logger.warning,
        ),
    )

# CLI entry point 
@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("config_path",
                type=click.Path(exists=True, dir_okay=False, path_type=Path))

def main(config_path: Path) -> None:
    """CLI entry: run the 3-step pipeline (lookup → load+annotate → export) from a YAML config.
    """
    # Configure logging
    logging.basicConfig(
        level=os.environ.get("MSPECVIZ_LOG_LEVEL", "INFO").upper(), # set MSPECVIZ_LOG_LEVEL=ERROR for a quiet run
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    start = time.time()
    config = _load_config(str(config_path))

    # Step 1: look up peptidoforms in idXMLs, bail if nothing matched
    grouped_psms, total_psms = _run_lookup(config)
    if total_psms == 0:
        sys.exit("Error: No PSMs found.")

    # Step 2: stream + annotate spectra against the lookup results
    psm_frame, psm_lookup = build_psm_table(grouped_psms)
    annotated_queries, reference_store, comparison_rows = _run_load(
        config, grouped_psms, psm_lookup, total_psms,
    )
    if not annotated_queries:
        sys.exit("Error: No spectra could be loaded or annotated.")

    # Step 3: export plots, tables, settings; report total wall time
    result = _run_export(
        config, annotated_queries, reference_store, grouped_psms,
        comparison_rows, psm_frame,
    )
    logger.info(f"Done in {time.time() - start:.1f}s. Exported to {result.export_root}")


if __name__ == "__main__":
    main()
