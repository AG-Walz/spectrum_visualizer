#!/usr/bin/env python3
"""
MSInsight CLI -- Terminal-based analysis pipeline.

Runs the full MSInsight analysis from the terminal without opening a browser.
All inputs are provided via a YAML config file (with optional CLI argument
overrides). Reuses the exact same backend functions as the Streamlit app.

Usage:
    python -m cli --config config.yaml
    python -m cli --config config.yaml --peptides peptides.txt --verbose
    python -m cli --config config.yaml --mz-tolerance 0.05 --no-plots --quiet

Exit codes:
    0 = success
    1 = pipeline failure
    2 = configuration / input error
"""
import argparse
import copy
import glob as globmod
import json
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from mspecviz.export_dirs import find_or_create_export_directory, commit_export_directory, sanitize_filename
from mspecviz.mirror_plot import (
    create_mirror_plot, create_single_spectrum_plot,
    figure_to_bytes, create_plotly_mirror_plot,
)
from mspecviz.insight_dataframes import (
    build_psm_table, build_peaks_dataframe,
    build_sequence_dataframe, build_heatmap_dataframe,
)
from mspecviz.peptide_search import search_peptides_in_idxmls, get_unique_peptides_from_idxmls
from mspecviz.query_processing import process_query_result, build_plot_tasks
from mspecviz.spectrum_loader import load_spectra_streaming
from mspecviz.summary_plots import generate_all_summary_plots
from mspecviz.common.text_utils import format_time, get_nested, parse_peptide_list

logger = logging.getLogger("msinsight.cli")


class TerminalProgress:
    """Simple terminal progress indicator with elapsed / ETA."""

    def __init__(self, label: str, total: int, quiet: bool = False):
        self.label = label
        self.total = max(total, 1)
        self.current = 0
        self.start = time.time()
        self.quiet = quiet

    def update(self, current: int, message: str = ""):
        self.current = current
        if self.quiet:
            return
        elapsed = time.time() - self.start
        pct = min(self.current / self.total * 100, 100)
        elapsed_str = format_time(elapsed)
        eta_str = ""
        if self.current > 0 and self.current < self.total:
            remaining = elapsed / self.current * (self.total - self.current)
            eta_str = f", ~{format_time(remaining)} remaining"
        extra = f"  {message}" if message else ""
        print(f"\r  {self.label} [{self.current}/{self.total}] "
              f"{pct:5.1f}%  ({elapsed_str} elapsed{eta_str}){extra}",
              end="", flush=True)

    def finish(self, message: str = ""):
        self.current = self.total
        if not self.quiet:
            elapsed = time.time() - self.start
            extra = f"  {message}" if message else ""
            print(f"\r  {self.label} [{self.total}/{self.total}] "
                  f"100.0%  ({format_time(elapsed)} elapsed){extra}")


@dataclass
class CLIConfig:
    """All pipeline parameters merged from settings.json, YAML, and CLI args."""

    # Input files (resolved absolute paths after glob expansion)
    idxml_files: List[str] = field(default_factory=list)
    mzml_files: List[str] = field(default_factory=list)

    # Peptides
    peptides: List[str] = field(default_factory=list)
    peptides_mode: str = ""  # "all", "file", or "inline"

    # Reference
    reference_substring: str = ""

    # Annotation
    ion_types: List[str] = field(default_factory=lambda: ["b", "y"])
    mz_tolerance: float = 0.02
    add_losses: bool = False
    add_first_prefix_ion: bool = True
    remove_precursor_peaks: bool = False
    intensity_mode: str = "Relative"

    # Export
    export_dir: str = "./exports"
    project_name: str = "MSInsight_analysis"
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    dpi: int = 300
    figsize: Tuple[int, int] = (12, 8)
    plot_unpaired: bool = False
    no_plots: bool = False

    # Performance
    workers: int = 4

    # Koina
    koina_enabled: bool = True
    koina_model: str = "ms2pip_timsTOF2024"
    koina_ce: int = 35
    koina_instrument: str = "TIMSTOF"


def load_settings_defaults() -> dict:
    """Load defaults from settings.json."""
    settings_path = PROJECT_ROOT / "settings.json"
    if settings_path.exists():
        with open(settings_path) as f:
            return json.load(f)
    return {}


def load_yaml_config(path: str) -> dict:
    """Load YAML config file."""
    config_path = Path(path)
    if not config_path.exists():
        print(f"Error: config file not found: {path}", file=sys.stderr)
        sys.exit(2)
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        print(f"Error: config file is not a valid YAML mapping: {path}", file=sys.stderr)
        sys.exit(2)
    return cfg


def expand_globs(patterns: List[str]) -> List[str]:
    """Expand glob patterns to sorted list of unique absolute paths."""
    result = []
    seen = set()
    for pattern in patterns:
        matches = sorted(globmod.glob(pattern, recursive=True))
        if not matches:
            logger.warning(f"No files matched pattern: {pattern}")
        for m in matches:
            abspath = str(Path(m).resolve())
            if abspath not in seen:
                seen.add(abspath)
                result.append(abspath)
    return result


def load_peptides(peptide_spec: str, idxml_files: List[str], workers: int,
                  quiet: bool) -> Tuple[List[str], str]:
    """
    Load peptides from spec (path, "all", or inline list in YAML).

    Returns (peptide_list, mode_string).
    """
    if isinstance(peptide_spec, list):
        # Inline list from YAML
        peps = [str(p).strip().upper() for p in peptide_spec if str(p).strip()]
        return peps, "inline"

    peptide_spec = str(peptide_spec).strip()

    if peptide_spec.lower() == "all":
        if not quiet:
            print("  Extracting all peptides from idXML files...")
        progress = TerminalProgress("Scanning", len(idxml_files), quiet)

        def on_progress(current, total, msg):
            progress.update(current, msg)

        peps = get_unique_peptides_from_idxmls(
            idxml_files, limit=10000,
            progress_callback=on_progress, workers=workers,
        )
        progress.finish(f"{len(peps)} unique peptides")
        return peps, "all"

    # Check if it is a file path
    pep_path = Path(peptide_spec)
    if pep_path.exists() and pep_path.is_file():
        text = pep_path.read_text()
        peps = parse_peptide_list(text)
        return peps, "file"

    # Last resort: treat as comma-separated list
    peps = [p.strip().upper() for p in peptide_spec.replace(",", "\n").split("\n") if p.strip()]
    return peps, "inline"


def build_config(yaml_cfg: dict, args: argparse.Namespace) -> CLIConfig:
    """Merge settings.json defaults, YAML config, and CLI overrides."""
    defaults = load_settings_defaults()
    cfg = CLIConfig()

    # --- Defaults from settings.json ---
    cfg.ion_types = get_nested(defaults, "ion_types", "default", default=["b", "y"])
    cfg.mz_tolerance = get_nested(defaults, "mz_tolerance", "default", default=0.02)
    cfg.add_losses = get_nested(defaults, "annotation", "add_losses", default=False)
    cfg.add_first_prefix_ion = get_nested(defaults, "annotation", "add_first_prefix_ion", default=True)
    cfg.remove_precursor_peaks = get_nested(defaults, "annotation", "remove_precursor_peaks", default=False)
    cfg.export_dir = get_nested(defaults, "export", "default_directory", default="./exports")
    cfg.project_name = get_nested(defaults, "export", "project_name", default="MSInsight_analysis")
    cfg.export_formats = get_nested(defaults, "export", "formats", default=["png"])
    cfg.dpi = get_nested(defaults, "export", "dpi", default=300)
    figsize = get_nested(defaults, "export", "figsize", default=[12, 8])
    cfg.figsize = tuple(figsize) if isinstance(figsize, list) else (12, 8)
    cfg.workers = get_nested(defaults, "performance", "workers", default=4)
    cfg.koina_enabled = get_nested(defaults, "koina", "enabled", default=True)
    cfg.koina_model = get_nested(defaults, "koina", "model", default="ms2pip_timsTOF2024")
    cfg.koina_ce = get_nested(defaults, "koina", "collision_energy", default=35)
    cfg.koina_instrument = get_nested(defaults, "koina", "instrument_type", default="TIMSTOF")

    # --- YAML overrides ---
    cfg.idxml_files = expand_globs(get_nested(yaml_cfg, "input", "idxml_files", default=[]))
    cfg.mzml_files = expand_globs(get_nested(yaml_cfg, "input", "mzml_files", default=[]))
    cfg.reference_substring = get_nested(yaml_cfg, "reference", "substring", default="")

    _yaml_mappings = [
        ("annotation", [
            ("ion_types", "ion_types"), ("mz_tolerance", "mz_tolerance", float),
            ("add_losses", "add_losses", bool), ("add_first_prefix_ion", "add_first_prefix_ion", bool),
            ("remove_precursor_peaks", "remove_precursor_peaks", bool),
            ("intensity_mode", "intensity_mode"),
        ]),
        ("export", [
            ("directory", "export_dir"), ("project_name", "project_name"),
            ("formats", "export_formats"), ("dpi", "dpi", int),
            ("figsize", "figsize", tuple), ("plot_unpaired", "plot_unpaired", bool),
        ]),
        ("performance", [("workers", "workers", int)]),
        ("koina", [
            ("enabled", "koina_enabled", bool), ("model", "koina_model"),
            ("collision_energy", "koina_ce", int), ("instrument_type", "koina_instrument"),
        ]),
    ]
    for section_name, keys in _yaml_mappings:
        section = yaml_cfg.get(section_name, {})
        for entry in keys:
            yaml_key, attr = entry[0], entry[1]
            cast = entry[2] if len(entry) > 2 else None
            if yaml_key in section:
                val = cast(section[yaml_key]) if cast else section[yaml_key]
                setattr(cfg, attr, val)

    # --- CLI overrides (highest priority) ---
    if args.idxml:
        cfg.idxml_files = expand_globs(args.idxml)
    if args.mzml:
        cfg.mzml_files = expand_globs(args.mzml)
    if args.reference_substring is not None:
        cfg.reference_substring = args.reference_substring
    if args.ion_types:
        cfg.ion_types = args.ion_types
    if args.mz_tolerance is not None:
        cfg.mz_tolerance = args.mz_tolerance
    if args.export_dir:
        cfg.export_dir = args.export_dir
    if args.project_name:
        cfg.project_name = args.project_name
    if args.formats:
        cfg.export_formats = args.formats
    if args.workers is not None:
        cfg.workers = args.workers
    if args.no_plots:
        cfg.no_plots = True
    if args.no_koina:
        cfg.koina_enabled = False
    if args.koina_model:
        cfg.koina_model = args.koina_model
    if args.koina_ce is not None:
        cfg.koina_ce = args.koina_ce
    if args.koina_instrument:
        cfg.koina_instrument = args.koina_instrument

    return cfg


def run_pipeline(cfg: CLIConfig, quiet: bool = False) -> int:
    """
    Execute the full MSInsight analysis pipeline.

    Returns exit code: 0 = success, 1 = failure.
    """
    pipeline_start = time.time()

    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    if not cfg.idxml_files:
        print("Error: No idXML files found. Check your config or --idxml argument.", file=sys.stderr)
        return 2
    if not cfg.mzml_files:
        print("Error: No mzML files found. Check your config or --mzml argument.", file=sys.stderr)
        return 2

    for label, files in [("idXML", cfg.idxml_files), ("mzML", cfg.mzml_files)]:
        for f in files:
            if not Path(f).exists():
                print(f"Error: {label} file not found: {f}", file=sys.stderr)
                return 2

    # ------------------------------------------------------------------
    # 2. Print configuration summary
    # ------------------------------------------------------------------
    if not quiet:
        print("\n" + "=" * 60)
        print("  MSInsight CLI Pipeline")
        print("=" * 60)
        print(f"  idXML files:    {len(cfg.idxml_files)}")
        for f in cfg.idxml_files:
            print(f"                  {Path(f).name}")
        print(f"  mzML files:     {len(cfg.mzml_files)}")
        for f in cfg.mzml_files:
            print(f"                  {Path(f).name}")
        print(f"  Peptides:       {len(cfg.peptides)} ({cfg.peptides_mode})")
        if len(cfg.peptides) <= 10:
            for p in cfg.peptides:
                print(f"                  {p}")
        else:
            for p in cfg.peptides[:5]:
                print(f"                  {p}")
            print(f"                  ... and {len(cfg.peptides) - 5} more")
        print(f"  Reference:      {cfg.reference_substring or '(none)'}")
        print(f"  Ion types:      {', '.join(cfg.ion_types)}")
        print(f"  m/z tolerance:  {cfg.mz_tolerance} Da")
        if cfg.koina_enabled:
            print(f"  Koina:          enabled ({cfg.koina_model})")
        else:
            print(f"  Koina:          disabled")
        print(f"  Workers:        {cfg.workers}")
        print(f"  Export dir:     {cfg.export_dir}")
        print(f"  Formats:        {', '.join(cfg.export_formats)}")
        print(f"  Plots:          {'skip' if cfg.no_plots else 'generate'}")
        print("=" * 60)
        print()

    # ------------------------------------------------------------------
    # 3. Search peptides in idXML files
    # ------------------------------------------------------------------
    if not quiet:
        print("Step 1: Searching peptides in idXML files...")

    search_progress = TerminalProgress("Searching", len(cfg.idxml_files), quiet)
    search_errors = []

    def search_progress_cb(current, total, msg):
        search_progress.update(current, msg)

    def search_error_cb(msg):
        search_errors.append(msg)
        logger.warning(msg)

    variant_groups = search_peptides_in_idxmls(
        idxml_paths=cfg.idxml_files,
        target_peptides=cfg.peptides,
        mzml_paths=cfg.mzml_files,
        progress_callback=search_progress_cb,
        error_callback=search_error_cb,
        workers=min(cfg.workers, len(cfg.idxml_files)),
    )
    search_progress.finish()

    if search_errors and not quiet:
        for err in search_errors:
            print(f"  Warning: {err}")

    total_psms = sum(len(variant_group.matches) for variant_group in variant_groups.values())
    if not quiet:
        print(f"  Found {total_psms} PSMs across {len(variant_groups)} variant(s)\n")

    if total_psms == 0:
        print("Error: No PSMs found. Check your peptide list and idXML files.", file=sys.stderr)
        return 1

    # Build PSM table and lookup
    vg_copy = copy.deepcopy(variant_groups)
    psm_lf, psm_lookup = build_psm_table(vg_copy)

    # ------------------------------------------------------------------
    # 4. Stream spectra from mzML files
    # ------------------------------------------------------------------
    if not quiet:
        print("Step 2: Loading and annotating spectra...")

    normalize = cfg.intensity_mode == "Relative"
    annotation_kwargs = dict(
        ion_types=cfg.ion_types,
        mz_tolerance=cfg.mz_tolerance,
        normalize=normalize,
        add_losses=cfg.add_losses,
        add_first_prefix_ion=cfg.add_first_prefix_ion,
        remove_precursor_peaks=cfg.remove_precursor_peaks,
    )

    annotated_queries = []        # (psm_id, LoadedSpectrum, annotated_df)
    reference_store = {}          # variant_key -> [(LoadedSpectrum, annotated_df)]
    comparison_rows = []          # dicts for comparison summary

    load_progress = TerminalProgress("Loading", total_psms, quiet)
    loaded_count = [0]

    def on_reference_loaded(variant_key, count):
        pass

    def on_query_ready(result):
        """Process each query: annotate, pair with references, compute metrics."""
        processed = process_query_result(
            result, psm_lookup, annotation_kwargs, reference_store,
            koina_enabled=cfg.koina_enabled,
            koina_model=cfg.koina_model,
            koina_ce=cfg.koina_ce,
            koina_instrument=cfg.koina_instrument,
        )
        if processed is None:
            return

        psm_id, query_spec, query_df, comp_rows = processed
        annotated_queries.append((psm_id, query_spec, query_df))
        comparison_rows.extend(comp_rows)

        loaded_count[0] += 1
        load_progress.update(loaded_count[0])

    def on_file_start(path, is_reference):
        kind = "reference" if is_reference else "query"
        logger.info(f"Processing {kind}: {Path(path).name}")

    def on_file_complete(path, loaded, total):
        logger.info(f"Completed {Path(path).name}: {loaded}/{total} spectra")

    def error_callback(error):
        logger.error(f"Streaming error: {error}")

    load_spectra_streaming(
        variant_groups=vg_copy,
        mzml_paths=[Path(m) for m in cfg.mzml_files],
        reference_substring=cfg.reference_substring,
        on_reference_loaded=on_reference_loaded,
        on_query_ready=on_query_ready,
        on_file_start=on_file_start,
        on_file_complete=on_file_complete,
        error_callback=error_callback,
        workers=cfg.workers,
    )
    load_progress.finish()

    if not annotated_queries:
        print("Error: No spectra could be loaded or annotated.", file=sys.stderr)
        return 1

    if not quiet:
        print(f"  Annotated {len(annotated_queries)} spectra, "
              f"{len(comparison_rows)} comparison(s)\n")

    # ------------------------------------------------------------------
    # 5. Build Polars DataFrames
    # ------------------------------------------------------------------
    if not quiet:
        print("Step 3: Building DataFrames...")

    peaks_lf = build_peaks_dataframe([(psm_id, df) for psm_id, _, df in annotated_queries])
    seq_lf = build_sequence_dataframe([(psm_id, s) for psm_id, s, _ in annotated_queries])
    heatmap_lf = build_heatmap_dataframe([s for _, s, _ in annotated_queries])

    comparison_df = pd.DataFrame(comparison_rows) if comparison_rows else pd.DataFrame()

    if not quiet:
        peaks_count = peaks_lf.select("x").collect().height
        print(f"  PSMs: {psm_lf.select('id').collect().height}, "
              f"Peaks: {peaks_count}, "
              f"Comparisons: {len(comparison_df)}\n")

    # ------------------------------------------------------------------
    # 6. Create export directory and save run settings
    # ------------------------------------------------------------------
    if not quiet:
        print("Step 4: Exporting results...")

    run_settings = {
        "idxml_files": cfg.idxml_files,
        "mzml_files": cfg.mzml_files,
        "peptides": cfg.peptides,
        "reference_substring": cfg.reference_substring,
        "ion_types": cfg.ion_types,
        "mz_tolerance": cfg.mz_tolerance,
        "add_losses": cfg.add_losses,
        "add_first_prefix_ion": cfg.add_first_prefix_ion,
        "remove_precursor_peaks": cfg.remove_precursor_peaks,
        "intensity_mode": cfg.intensity_mode,
        "koina_enabled": cfg.koina_enabled,
        "koina_model": cfg.koina_model,
        "koina_ce": cfg.koina_ce,
        "koina_instrument": cfg.koina_instrument,
        "export_formats": cfg.export_formats,
        "dpi": cfg.dpi,
        "figsize": list(cfg.figsize),
    }

    root_dir, png_dir, svg_dir, pdf_dir = find_or_create_export_directory(
        cfg.export_dir, cfg.project_name, run_settings,
    )
    commit_export_directory(root_dir, run_settings)

    # Create format subdirectories as needed
    for fmt in cfg.export_formats:
        fmt_dir = {"png": png_dir, "svg": svg_dir, "pdf": pdf_dir}.get(fmt)
        if fmt_dir:
            fmt_dir.mkdir(parents=True, exist_ok=True)

    saved_files = []

    # ------------------------------------------------------------------
    # 7. Export Parquet DataFrames
    # ------------------------------------------------------------------
    for name, lf in [("psm", psm_lf), ("peaks", peaks_lf),
                     ("seq", seq_lf), ("heatmap", heatmap_lf)]:
        out = root_dir / f"{name}.parquet"
        lf.collect().write_parquet(out)
        saved_files.append(str(out))

    if not quiet:
        print(f"  Saved 4 Parquet files")

    # ------------------------------------------------------------------
    # 8. Export comparison summary CSV/XLSX
    # ------------------------------------------------------------------
    if not comparison_df.empty:
        csv_path = root_dir / "comparison_summary.csv"
        comparison_df.to_csv(csv_path, index=False)
        saved_files.append(str(csv_path))

        try:
            xlsx_path = root_dir / "comparison_summary.xlsx"
            comparison_df.to_excel(xlsx_path, index=False, engine="openpyxl")
            saved_files.append(str(xlsx_path))
        except ImportError:
            logger.warning("openpyxl not installed -- skipping XLSX export")
        except Exception as e:
            logger.warning(f"XLSX export failed: {e}")

        if not quiet:
            print(f"  Saved comparison summary (CSV"
                  f"{'/XLSX' if (root_dir / 'comparison_summary.xlsx').exists() else ''})")

    # ------------------------------------------------------------------
    # 9. Generate mirror plots (parallel)
    # ------------------------------------------------------------------
    if not cfg.no_plots:
        n_plots = _export_mirror_plots(
            annotated_queries, reference_store, vg_copy,
            cfg, root_dir, png_dir, svg_dir, pdf_dir,
            normalize, quiet,
        )
        saved_files.extend(["(plots)"] if n_plots > 0 else [])

        if not quiet:
            print(f"  Generated {n_plots} plot file(s)")

        # Summary plots
        if not comparison_df.empty:
            summary_formats = [f for f in cfg.export_formats if f in ("png", "svg", "pdf")]
            if not summary_formats:
                summary_formats = ["png"]
            results = generate_all_summary_plots(
                source=comparison_df,
                output_dir=str(root_dir),
                formats=summary_formats,
            )
            n_summary = sum(len(v) for v in results.values())
            if not quiet:
                print(f"  Generated {n_summary} summary plot(s)")

    # ------------------------------------------------------------------
    # 10. Completion summary
    # ------------------------------------------------------------------
    elapsed = time.time() - pipeline_start
    if not quiet:
        print()
        print("=" * 60)
        print("  Pipeline Complete")
        print("=" * 60)
        print(f"  Variants:     {len(variant_groups)}")
        print(f"  PSMs:         {total_psms}")
        print(f"  Spectra:      {len(annotated_queries)}")
        print(f"  Comparisons:  {len(comparison_rows)}")
        print(f"  Export dir:   {root_dir}")
        print(f"  Elapsed:      {format_time(elapsed)}")
        print("=" * 60)
        print()

    return 0


def _generate_single_plot(task: dict) -> int:
    """Worker function for parallel plot generation. Returns file count."""
    count = 0

    if task["type"] == "mirror":
        fig = create_mirror_plot(
            upper_spectrum=task["query_df"],
            lower_spectrum=task["ref_df"],
            title=task["title"],
            upper_label=task["query_label"],
            lower_label=task["ref_label"],
            metrics_text=task["metrics_text"],
            normalize=task["normalize"],
            figsize=task.get("figsize", (12, 8)),
            dpi=task.get("dpi", 300),
        )
    else:
        fig = create_single_spectrum_plot(
            spectrum=task["query_df"],
            title=task["title"],
            label=task.get("query_label", ""),
            normalize=task["normalize"],
            figsize=task.get("figsize", (12, 8)),
            dpi=task.get("dpi", 300),
        )

    for fmt_dir_key, fmt in [("png_dir", "png"), ("svg_dir", "svg"), ("pdf_dir", "pdf")]:
        fmt_dir = task.get(fmt_dir_key, "")
        if not fmt_dir:
            continue
        out_dir = Path(fmt_dir) / task["folder_name"]
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"{task['safe_name']}.{fmt}"
        data = figure_to_bytes(fig, format=fmt)
        out_file.write_bytes(data)
        count += 1

    plt.close(fig)

    # HTML interactive plot
    if task.get("html_dir") and task["type"] == "mirror":
        plotly_fig = create_plotly_mirror_plot(
            upper_spectrum=task["query_df"],
            lower_spectrum=task["ref_df"],
            title=task["title"],
            upper_label=task["query_label"],
            lower_label=task["ref_label"],
            metrics_text=task["metrics_text"],
            normalize=task["normalize"],
        )
        html_out_dir = Path(task["html_dir"]) / task["folder_name"]
        html_out_dir.mkdir(parents=True, exist_ok=True)
        html_path = html_out_dir / f"{task['safe_name']}.html"
        plotly_fig.write_html(str(html_path), include_plotlyjs="cdn")
        count += 1

    return count


def _export_mirror_plots(annotated_queries, reference_store, variant_groups,
                         cfg, root_dir, png_dir, svg_dir, pdf_dir,
                         normalize, quiet):
    """Generate mirror / single-spectrum plots using parallel workers."""
    # Determine which format directories are active
    fmt_dirs = {}
    for fmt, d in [("png", png_dir), ("svg", svg_dir), ("pdf", pdf_dir)]:
        if fmt in cfg.export_formats:
            fmt_dirs[fmt] = str(d)

    html_dir = str(root_dir / "html") if "html" in cfg.export_formats else ""

    tasks = build_plot_tasks(
        annotated_queries, reference_store, variant_groups,
        normalize=normalize,
        sanitize_fn=sanitize_filename,
        fmt_dirs=fmt_dirs,
        figsize=cfg.figsize,
        dpi=cfg.dpi,
        html_dir=html_dir,
        plot_unpaired=cfg.plot_unpaired,
    )

    if not tasks:
        return 0

    plot_progress = TerminalProgress("Generating plots", len(tasks), quiet)
    total_files = 0
    effective_workers = min(cfg.workers, len(tasks))

    if effective_workers <= 1:
        for i, task in enumerate(tasks):
            total_files += _generate_single_plot(task)
            plot_progress.update(i + 1)
    else:
        completed = 0
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {executor.submit(_generate_single_plot, t): t for t in tasks}
            for future in as_completed(futures):
                completed += 1
                try:
                    total_files += future.result()
                except Exception as e:
                    logger.error(f"Plot generation failed: {e}")
                plot_progress.update(completed)

    plot_progress.finish()
    return total_files


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="msinsight-cli",
        description="MSInsight CLI -- Terminal-based MS/MS spectrum analysis pipeline.",
        epilog="Example: python -m cli --config cli/config/config_example.yaml --verbose",
    )

    # Required
    parser.add_argument("-c", "--config", required=True,
                        help="Path to YAML config file")

    # Input overrides
    parser.add_argument("--idxml", nargs="+", metavar="PATTERN",
                        help="idXML file paths or glob patterns (overrides config)")
    parser.add_argument("--mzml", nargs="+", metavar="PATTERN",
                        help="mzML file paths or glob patterns (overrides config)")
    parser.add_argument("--peptides", metavar="PATH_OR_ALL",
                        help='Path to peptides.txt or "all" (overrides config)')
    parser.add_argument("--reference-substring", metavar="STR",
                        help="Reference file substring (overrides config)")

    # Annotation overrides
    parser.add_argument("--ion-types", nargs="+", metavar="ION",
                        help="Ion types, e.g. b y (overrides config)")
    parser.add_argument("--mz-tolerance", type=float, metavar="DA",
                        help="m/z tolerance in Da (overrides config)")

    # Export overrides
    parser.add_argument("--export-dir", metavar="DIR",
                        help="Export directory (overrides config)")
    parser.add_argument("--project-name", metavar="NAME",
                        help="Project name (overrides config)")
    parser.add_argument("--formats", nargs="+", metavar="FMT",
                        choices=["png", "svg", "pdf", "html"],
                        help="Export formats: png svg pdf html (overrides config)")
    parser.add_argument("--workers", type=int, metavar="N",
                        help="Parallel workers (overrides config)")

    # Flags
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plot generation, export CSV/Parquet only")
    parser.add_argument("--no-koina", action="store_true",
                        help="Disable Koina AI fallback predictions")

    # Koina overrides
    parser.add_argument("--koina-model", metavar="ID",
                        help="Koina model ID (overrides config)")
    parser.add_argument("--koina-ce", type=int, metavar="EV",
                        help="Koina collision energy (overrides config)")
    parser.add_argument("--koina-instrument", metavar="TYPE",
                        help="Koina instrument type (overrides config)")

    # Verbosity
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output (DEBUG logging)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress non-error output")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Configure logging
    if args.quiet:
        log_level = logging.ERROR
    elif args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    # Load YAML config
    yaml_cfg = load_yaml_config(args.config)

    # Build merged config
    cfg = build_config(yaml_cfg, args)

    # Load peptides (needs to happen after config is built, before pipeline)
    peptide_spec = args.peptides if args.peptides else yaml_cfg.get("peptides", [])
    if peptide_spec:
        cfg.peptides, cfg.peptides_mode = load_peptides(
            peptide_spec, cfg.idxml_files, cfg.workers, args.quiet,
        )
    else:
        print("Error: No peptides specified. Use 'peptides' in config or --peptides.", file=sys.stderr)
        sys.exit(2)

    if not cfg.peptides:
        print("Error: Peptide list is empty.", file=sys.stderr)
        sys.exit(2)

    # Run pipeline
    exit_code = run_pipeline(cfg, quiet=args.quiet)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
