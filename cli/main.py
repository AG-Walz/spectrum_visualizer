"""MSpecViz CLI entry: run the 3-step pipeline (search → load+annotate → export) from a YAML config; 
invoke via ``python -m cli <config.yaml>``."""
import copy
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

from core.common.defaults import DEFAULTS  # noqa: E402
from core.dataframes import build_psm_table, build_visualization_frames  # noqa: E402
from core.export import (  # noqa: E402
    ExportCallbacks, ExportSelection, PlotSettings,
    build_run_settings, run_export,
)
from core.pipeline import (  # noqa: E402
    build_annotation_kwargs, build_koina_kwargs, run_streaming_annotation,
)
from core.search import search_peptides_in_idxmls  # noqa: E402

logger = logging.getLogger("core.cli")

# CLI-only fields not in DEFAULTS — populated by the user's YAML.
_CLI_EXTRAS: dict = {
    "idxml_files": [], "mzml_files": [],
    "peptide_source": "", "peptides": [],
    "reference_substring": "", "no_plots": False,
}


def _expand_globs(patterns: list) -> list:
    """Expand shell-glob ``patterns`` (e.g. ``"runs/**/ProT*.idXML"``) into a sorted, unique list of absolute paths."""
    result, seen = [], set()
    for pattern in patterns:
        for match in sorted(globmod.glob(pattern, recursive=True)):
            absolute = str(Path(match).resolve())
            if absolute not in seen:
                seen.add(absolute)
                result.append(absolute)
    return result


def _load_peptides(source) -> list:
    """Resolve ``source`` (YAML list, file path, or raw comma/newline string) into a list of upper-cased sequences."""
    # Already a list: strip and uppercase each entry
    if isinstance(source, list):
        return [str(peptide).strip().upper() for peptide in source if str(peptide).strip()]
    # File path: read its contents; otherwise treat the string itself as the peptide blob
    source = str(source).strip()
    text = Path(source).read_text() if Path(source).is_file() else source
    return [line.strip().upper() for line in text.replace(",", "\n").splitlines() if line.strip()]


def _load_config(path: str) -> SimpleNamespace:
    """Parse the YAML at ``path``, overlay on :data:`DEFAULTS`, expand file globs, and resolve the peptide source into a ``SimpleNamespace``."""
    # Parse: read YAML file, error out on non-mapping payload
    config_path = Path(path)
    yaml_config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(yaml_config, dict):
        sys.exit(f"Error: config must be a YAML mapping: {path}")

    # Merge: CLI extras + DEFAULTS overlaid by user YAML (last write wins)
    merged: dict = {**_CLI_EXTRAS, **asdict(DEFAULTS), **yaml_config}
    config = SimpleNamespace(**merged)

    # Expand input globs and require non-empty idXML and mzML lists
    config.idxml_files = _expand_globs(config.idxml_files)
    config.mzml_files = _expand_globs(config.mzml_files)
    if not config.idxml_files or not config.mzml_files:
        sys.exit("Error: idxml_files and mzml_files must both be set.")
    if not config.peptide_source:
        sys.exit("Error: 'peptide_source' is required.")
    config.peptides = _load_peptides(config.peptide_source)
    if not config.peptides:
        sys.exit("Error: peptide list is empty.")
    return config


def _run_search(config) -> tuple:
    """Step 1: scan idXMLs for target peptides with a tqdm bar; return ``(variant_groups, total_psms)``."""
    logger.info("Step 1: Searching peptides in idXML files")
    bar = tqdm(total=len(config.idxml_files), desc="Searching", unit="file")
    errors: list = []
    variant_groups = search_peptides_in_idxmls(
        idxml_paths=config.idxml_files,
        target_peptides=config.peptides,
        mzml_paths=config.mzml_files,
        progress_callback=lambda current, total, message: bar.update(current - bar.n),
        error_callback=errors.append,
        workers=min(config.workers, len(config.idxml_files)),
    )
    bar.close()
    for error in errors:
        logger.warning(error)
    total_psms = sum(len(psm_list) for psm_list in variant_groups.values())
    logger.info(f"Found {total_psms} PSMs across {len(variant_groups)} variants")
    return variant_groups, total_psms


def _run_load(config, variant_groups, psm_lookup, total_psms) -> tuple:
    """Step 2: stream spectra from mzMLs and annotate fragment ions; return ``(annotated_queries, reference_store, comparison_rows)``."""
    logger.info("Step 2: Loading and annotating spectra")
    getter = lambda key, default: getattr(config, key, default)
    bar = tqdm(total=total_psms, desc="Loading", unit="spec")
    annotated_queries, reference_store, comparison_rows = run_streaming_annotation(
        variant_groups=variant_groups,
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


def _run_export(config, annotated_queries, reference_store, variant_groups,
                comparison_rows, psm_frame):
    """Step 3: write parquet/CSV tables and (unless ``--no-plots``) mirror + summary plots via :func:`core.export.run_export`."""
    logger.info("Step 3: Exporting results")
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
        variant_groups=variant_groups,
        comparison_frame=pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame(),
        selection=ExportSelection(
            parquet=True, csv=True,
            mirror_plots=not config.no_plots,
            summary_plots=not config.no_plots,
            html="html" in config.export_formats,
        ),
        plot_settings=PlotSettings(
            export_formats=list(config.export_formats),
            figure_size=tuple(config.figure_size),
            export_resolution=config.export_resolution,
            plot_unpaired_queries=config.plot_unpaired_queries,
            workers=config.workers,
        ),
        callbacks=ExportCallbacks(
            on_step=logger.info, on_plot_progress=_on_plot_progress,
            on_warning=logger.warning,
        ),
    )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("config_path",
                type=click.Path(exists=True, dir_okay=False, path_type=Path))
def main(config_path: Path) -> None:
    """CLI entry: run the 3-step pipeline (search → load+annotate → export) from a YAML config; set ``MSPECVIZ_LOG_LEVEL=ERROR`` for a quiet run."""
    logging.basicConfig(
        level=os.environ.get("MSPECVIZ_LOG_LEVEL", "INFO").upper(),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    start = time.time()
    config = _load_config(str(config_path))

    # Step 1: search idXMLs for peptide variants, bail if nothing matched
    variant_groups, total_psms = _run_search(config)
    if total_psms == 0:
        sys.exit("Error: No PSMs found.")

    # Step 2: stream + annotate spectra against the search results
    variant_groups_copy = copy.deepcopy(variant_groups)
    psm_frame, psm_lookup = build_psm_table(variant_groups_copy)
    annotated_queries, reference_store, comparison_rows = _run_load(
        config, variant_groups_copy, psm_lookup, total_psms,
    )
    if not annotated_queries:
        sys.exit("Error: No spectra could be loaded or annotated.")

    # Step 3: export plots, tables, settings; report total wall time
    result = _run_export(
        config, annotated_queries, reference_store, variant_groups_copy,
        comparison_rows, psm_frame,
    )
    logger.info(f"Done in {time.time() - start:.1f}s. Exported to {result.export_root}")


if __name__ == "__main__":
    main()
