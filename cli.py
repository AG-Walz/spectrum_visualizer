#!/usr/bin/env python3
"""
MSpecViz CLI — Terminal-based analysis pipeline.

Runs the full MSpecViz analysis from the terminal without opening a browser.
All inputs are provided via a YAML config file (with optional CLI argument
overrides). Reuses the exact same backend functions as the Streamlit app.

Usage:
    python cli.py --config config.yaml
    python cli.py --config config.yaml --peptides peptides.txt --verbose
    python cli.py --config config.yaml --mz-tolerance 0.05 --add-losses --quiet
"""
import argparse
import copy
import glob as globmod
import json
import logging
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

# ---------------------------------------------------------------------------
# Path setup — ensure project root is importable
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logger = logging.getLogger("mspecviz.cli")


# ============================================================================
# CLIConfig — dataclass holding all pipeline parameters
# ============================================================================

@dataclass
class CLIConfig:
    """All pipeline parameters merged from settings.json, YAML, and CLI args."""

    # Input files (resolved, absolute paths)
    idxml_files: List[str] = field(default_factory=list)
    mzml_files: List[str] = field(default_factory=list)
    parquet_file: Optional[str] = None

    # Peptides
    peptides: List[str] = field(default_factory=list)

    # Reference
    reference_substring: str = ""

    # Annotation
    ion_types: List[str] = field(default_factory=lambda: ["b", "y"])
    mz_tolerance: float = 0.02
    add_losses: bool = False
    add_first_prefix_ion: bool = True
    remove_precursor_peaks: bool = False
    intensity_mode: str = "Relative"  # "Relative" or "Raw"

    # Export
    export_dir: str = "./exports"
    project_name: str = "MSpecViz_analysis"
    export_formats: List[str] = field(default_factory=lambda: ["png"])
    dpi: int = 300
    figsize: Tuple[int, int] = (12, 8)
    plot_unpaired: bool = False

    # Cache
    cache_enabled: bool = True
    cache_subdirectory: str = ".mspecviz_cache"

    # Parallelism
    workers: int = 4

    # tmpfs staging
    tmpfs_staging: bool = True
    tmpfs_path: str = "/dev/shm"
    tmpfs_ram_fraction: float = 0.40

    # Verbosity
    verbose: bool = False
    quiet: bool = False

    @property
    def normalize(self) -> bool:
        return self.intensity_mode == "Relative"


def _expand_file_globs(patterns: List[str]) -> List[str]:
    """Expand glob patterns into absolute file paths, preserving order.

    Supports recursive patterns with ** (e.g. /path/**/*.idXML).
    """
    result = []
    seen = set()
    for pattern in patterns:
        recursive = "**" in pattern
        expanded = sorted(globmod.glob(pattern, recursive=recursive))
        if not expanded:
            # Not a glob, treat as literal path
            p = str(Path(pattern).resolve())
            if p not in seen:
                result.append(p)
                seen.add(p)
        else:
            for f in expanded:
                p = str(Path(f).resolve())
                if p not in seen:
                    result.append(p)
                    seen.add(p)
    return result


def _strip_modifications(seq: str) -> str:
    """Strip parenthetical modification annotations to get base sequence.

    E.g. 'ATYSVGM(Oxidation)QK' -> 'ATYSVGMQK'
    """
    import re
    return re.sub(r'\([^)]*\)', '', seq)


def _load_peptides(source) -> List[str]:
    """
    Load peptides from a .txt file path, inline list, or the keyword "all".

    Accepts:
      - "all" — special keyword, returns empty list (resolved later from idXMLs)
      - A string path to a .txt file (one peptide per line, no quotes)
      - A list of strings (from YAML inline list, for backwards compatibility)

    Modifications like (Oxidation) are stripped to get base sequences.
    Duplicate base sequences are removed.

    Returns list of uppercase, stripped base sequences.
    An empty list signals "all peptides" mode.
    """
    raw = []
    if isinstance(source, str):
        if source.strip().lower() == "all":
            return []  # sentinel: resolved later from idXMLs
        # Treat as file path
        p = Path(source)
        if not p.exists():
            raise FileNotFoundError(f"Peptide file not found: {source}")
        raw = p.read_text().strip().splitlines()
    elif isinstance(source, list):
        # Could be list of strings or list containing a file path
        if len(source) == 1 and Path(source[0]).suffix in ('.txt', '.csv', '.tsv'):
            return _load_peptides(source[0])
        raw = [str(p) for p in source]
    else:
        raise ValueError(f"Invalid peptides source: {type(source)}")

    # Strip modifications, uppercase, deduplicate while preserving order
    seen = set()
    result = []
    for line in raw:
        seq = _strip_modifications(line.strip()).upper()
        if seq and seq not in seen:
            result.append(seq)
            seen.add(seq)
    return result


def load_settings_json() -> dict:
    """Load settings.json from the project root."""
    settings_path = PROJECT_ROOT / "settings.json"
    if settings_path.exists():
        with open(settings_path) as f:
            return json.load(f)
    return {}


def _get_nested(d: dict, *keys, default=None):
    """Safely get nested dict value."""
    val = d
    for k in keys:
        if isinstance(val, dict):
            val = val.get(k)
        else:
            return default
    return val if val is not None else default


def load_config(yaml_path: str, cli_args: argparse.Namespace) -> CLIConfig:
    """
    Build CLIConfig by merging three layers (lowest to highest priority):
      1. settings.json defaults
      2. YAML config file values
      3. CLI argument overrides
    """
    # Layer 1: settings.json
    settings = load_settings_json()
    cfg = CLIConfig(
        ion_types=_get_nested(settings, "ion_types", "default", default=["b", "y"]),
        mz_tolerance=_get_nested(settings, "mz_tolerance", "default", default=0.02),
        add_losses=_get_nested(settings, "annotation", "add_losses", default=False),
        add_first_prefix_ion=_get_nested(settings, "annotation", "add_first_prefix_ion", default=True),
        remove_precursor_peaks=_get_nested(settings, "annotation", "remove_precursor_peaks", default=False),
        export_dir=_get_nested(settings, "export", "default_directory", default="./exports"),
        project_name=_get_nested(settings, "export", "project_name", default="MSpecViz_analysis"),
        export_formats=_get_nested(settings, "export", "formats", default=["png"]),
        dpi=_get_nested(settings, "export", "dpi", default=300),
        figsize=tuple(_get_nested(settings, "export", "figsize", default=[12, 8])),
        cache_enabled=_get_nested(settings, "cache", "enabled", default=True),
        cache_subdirectory=_get_nested(settings, "cache", "subdirectory", default=".mspecviz_cache"),
        workers=_get_nested(settings, "performance", "workers", default=4),
        tmpfs_staging=_get_nested(settings, "performance", "tmpfs_staging", default=True),
        tmpfs_path=_get_nested(settings, "performance", "tmpfs_path", default="/dev/shm"),
        tmpfs_ram_fraction=_get_nested(settings, "performance", "tmpfs_ram_fraction", default=0.40),
    )

    # Layer 2: YAML config file
    with open(yaml_path) as f:
        yaml_cfg = yaml.safe_load(f) or {}

    inp = yaml_cfg.get("input", {})
    if inp.get("idxml_files"):
        cfg.idxml_files = _expand_file_globs(inp["idxml_files"])
    if inp.get("mzml_files"):
        cfg.mzml_files = _expand_file_globs(inp["mzml_files"])
    if inp.get("parquet_file"):
        p = Path(inp["parquet_file"]).resolve()
        cfg.parquet_file = str(p) if p.exists() else None

    if yaml_cfg.get("peptides"):
        cfg.peptides = _load_peptides(yaml_cfg["peptides"])

    ref = yaml_cfg.get("reference", {})
    if "substring" in ref:
        cfg.reference_substring = ref["substring"] or ""

    ann = yaml_cfg.get("annotation", {})
    if "ion_types" in ann:
        cfg.ion_types = ann["ion_types"]
    if "mz_tolerance" in ann:
        cfg.mz_tolerance = float(ann["mz_tolerance"])
    if "add_losses" in ann:
        cfg.add_losses = bool(ann["add_losses"])
    if "add_first_prefix_ion" in ann:
        cfg.add_first_prefix_ion = bool(ann["add_first_prefix_ion"])
    if "remove_precursor_peaks" in ann:
        cfg.remove_precursor_peaks = bool(ann["remove_precursor_peaks"])
    if "intensity_mode" in ann:
        cfg.intensity_mode = ann["intensity_mode"]

    exp = yaml_cfg.get("export", {})
    if "directory" in exp:
        cfg.export_dir = exp["directory"]
    if "project_name" in exp:
        cfg.project_name = exp["project_name"]
    if "formats" in exp:
        cfg.export_formats = exp["formats"]
    if "dpi" in exp:
        cfg.dpi = int(exp["dpi"])
    if "figsize" in exp:
        cfg.figsize = tuple(exp["figsize"])
    if "plot_unpaired" in exp:
        cfg.plot_unpaired = bool(exp["plot_unpaired"])

    cache = yaml_cfg.get("cache", {})
    if "enabled" in cache:
        cfg.cache_enabled = bool(cache["enabled"])
    if "subdirectory" in cache:
        cfg.cache_subdirectory = cache["subdirectory"]

    perf = yaml_cfg.get("performance", {})
    if "workers" in perf:
        cfg.workers = int(perf["workers"])
    if "tmpfs_staging" in perf:
        cfg.tmpfs_staging = bool(perf["tmpfs_staging"])
    if "tmpfs_path" in perf:
        cfg.tmpfs_path = str(perf["tmpfs_path"])
    if "tmpfs_ram_fraction" in perf:
        cfg.tmpfs_ram_fraction = float(perf["tmpfs_ram_fraction"])

    # Layer 3: CLI argument overrides
    if cli_args.idxml:
        cfg.idxml_files = _expand_file_globs(cli_args.idxml)
    if cli_args.mzml:
        cfg.mzml_files = _expand_file_globs(cli_args.mzml)
    if cli_args.parquet:
        p = Path(cli_args.parquet).resolve()
        cfg.parquet_file = str(p)
    if cli_args.peptides:
        cfg.peptides = _load_peptides(cli_args.peptides)
    if cli_args.reference_substring is not None:
        cfg.reference_substring = cli_args.reference_substring
    if cli_args.ion_types:
        cfg.ion_types = cli_args.ion_types
    if cli_args.mz_tolerance is not None:
        cfg.mz_tolerance = cli_args.mz_tolerance
    if cli_args.add_losses:
        cfg.add_losses = True
    if cli_args.no_first_prefix_ion:
        cfg.add_first_prefix_ion = False
    if cli_args.remove_precursor_peaks:
        cfg.remove_precursor_peaks = True
    if cli_args.raw_intensity:
        cfg.intensity_mode = "Raw"
    if cli_args.export_dir:
        cfg.export_dir = cli_args.export_dir
    if cli_args.project_name:
        cfg.project_name = cli_args.project_name
    if cli_args.formats:
        cfg.export_formats = cli_args.formats
    if cli_args.dpi is not None:
        cfg.dpi = cli_args.dpi
    if cli_args.plot_unpaired:
        cfg.plot_unpaired = True
    if cli_args.no_cache:
        cfg.cache_enabled = False
    if cli_args.workers is not None:
        cfg.workers = cli_args.workers
    if cli_args.no_tmpfs:
        cfg.tmpfs_staging = False
    if cli_args.verbose:
        cfg.verbose = True
    if cli_args.quiet:
        cfg.quiet = True

    return cfg


# ============================================================================
# TerminalProgress — replaces Streamlit progress bars
# ============================================================================

class TerminalProgress:
    """Simple terminal progress indicator."""

    def __init__(self, total: int, label: str = "", quiet: bool = False):
        self.total = total
        self.label = label
        self.current = 0
        self.start_time = time.time()
        self.quiet = quiet

    def update(self, current: int = None):
        if self.quiet:
            return
        if current is not None:
            self.current = current
        else:
            self.current += 1

        elapsed = time.time() - self.start_time
        if self.total > 0:
            frac = min(self.current / self.total, 1.0)
            pct = frac * 100
            if frac > 0:
                remaining = (elapsed / frac) - elapsed
                remaining_str = format_time(remaining)
            else:
                remaining_str = "?"
            sys.stdout.write(
                f"\r  {self.label} [{self.current}/{self.total}] "
                f"{pct:5.1f}%  ({format_time(elapsed)} elapsed, ~{remaining_str} remaining)"
            )
        else:
            sys.stdout.write(f"\r  {self.label} {self.current} ({format_time(elapsed)} elapsed)")
        sys.stdout.flush()

    def finish(self, message: str = ""):
        if self.quiet:
            return
        elapsed = time.time() - self.start_time
        if message:
            sys.stdout.write(f"\r  {message} ({format_time(elapsed)})\n")
        else:
            sys.stdout.write(
                f"\r  {self.label} [{self.current}/{self.total}] Done ({format_time(elapsed)})\n"
            )
        sys.stdout.flush()


def format_time(seconds: float) -> str:
    """Format seconds as human-readable time."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def log_info(msg: str, quiet: bool = False):
    """Print info message unless quiet."""
    if not quiet:
        print(msg)


def log_error(msg: str):
    """Print error message to stderr."""
    print(f"ERROR: {msg}", file=sys.stderr)


def log_debug(msg: str, verbose: bool = False):
    """Print debug message if verbose."""
    if verbose:
        print(f"  DEBUG: {msg}")


# ============================================================================
# run_pipeline — core processing (translation of run_streaming_processing)
# ============================================================================

def run_pipeline(cfg: CLIConfig) -> int:
    """
    Run the full MSpecViz analysis pipeline.

    Returns exit code: 0 = success, 1 = processing failure.
    """
    from src.spectrum_loader import load_spectra_streaming, load_spectra_from_parquet, StreamingResult
    from src.peak_annotation import prepare_spectrum_for_plot
    from src.export import sanitize_filename, create_export_directory
    from src.spectrum_metrics import shorten_filenames

    # =========================================================================
    # Phase 0: Validate inputs
    # =========================================================================
    use_parquet = bool(cfg.parquet_file)

    if not use_parquet:
        if not cfg.idxml_files:
            log_error("No idXML files specified")
            return 1
        if not cfg.mzml_files:
            log_error("No mzML files specified")
            return 1

        # Validate files exist
        missing = [f for f in cfg.idxml_files if not Path(f).exists()]
        if missing:
            log_error(f"idXML files not found: {missing[:3]}{'...' if len(missing) > 3 else ''}")
            return 1
        missing = [f for f in cfg.mzml_files if not Path(f).exists()]
        if missing:
            log_error(f"mzML files not found: {missing[:3]}{'...' if len(missing) > 3 else ''}")
            return 1

        # "all" peptides mode: extract every unique peptide from idXMLs
        if not cfg.peptides:
            from src.peptide_search import get_unique_peptides_from_idxmls
            log_info("Extracting all peptides from idXML files...", cfg.quiet)
            cfg.peptides = get_unique_peptides_from_idxmls(cfg.idxml_files, limit=999999, workers=cfg.workers)
            if not cfg.peptides:
                log_error("No peptides found in idXML files")
                return 1
            log_info(f"Found {len(cfg.peptides)} unique peptides", cfg.quiet)
    else:
        if not Path(cfg.parquet_file).exists():
            log_error(f"Parquet file not found: {cfg.parquet_file}")
            return 1

    # =========================================================================
    # Print configuration summary
    # =========================================================================
    log_info("=" * 60, cfg.quiet)
    log_info("MSpecViz CLI Analysis", cfg.quiet)
    log_info("=" * 60, cfg.quiet)

    if use_parquet:
        log_info(f"Input:      Parquet ({Path(cfg.parquet_file).name})", cfg.quiet)
    else:
        log_info(f"Input:      {len(cfg.idxml_files)} idXML, {len(cfg.mzml_files)} mzML files", cfg.quiet)
        log_info(f"Peptides:   {', '.join(cfg.peptides[:5])}{'...' if len(cfg.peptides) > 5 else ''}", cfg.quiet)

    if cfg.reference_substring:
        log_info(f"Reference:  substring='{cfg.reference_substring}'", cfg.quiet)
    else:
        log_info("Reference:  (none — single spectrum mode)", cfg.quiet)

    log_info(f"Ions:       {cfg.ion_types}, tolerance={cfg.mz_tolerance} Da", cfg.quiet)
    log_info(f"Export:     {cfg.export_dir} [{', '.join(cfg.export_formats)}]", cfg.quiet)
    log_info(f"Project:    {cfg.project_name}", cfg.quiet)
    log_info(f"Workers:    {cfg.workers}", cfg.quiet)
    log_info("-" * 60, cfg.quiet)

    # =========================================================================
    # Phase 1: Auto-cache detection
    # =========================================================================
    loaded_from_cache = False
    cache_dir = None
    cache_fingerprint = None

    partial_cache_path = None  # set if a partial (incomplete) cache is found

    if not use_parquet and cfg.cache_enabled:
        from src.parquet_cache import compute_cache_fingerprint, find_existing_cache

        cache_fingerprint = compute_cache_fingerprint(
            cfg.idxml_files, cfg.mzml_files, cfg.peptides, cfg.reference_substring
        )
        cache_dir = str(Path(cfg.export_dir) / cfg.cache_subdirectory)

        # Check for complete cache first
        cached_path = find_existing_cache(cache_dir, cache_fingerprint)
        if cached_path:
            cfg.parquet_file = cached_path
            use_parquet = True
            loaded_from_cache = True
            log_info(f"Cache hit! Loading from: {Path(cached_path).name}", cfg.quiet)
        else:
            # Check for partial (incomplete) cache
            partial_path = find_existing_cache(cache_dir, cache_fingerprint, allow_partial=True)
            if partial_path:
                from src.parquet_cache import get_parquet_row_counts
                p_refs, p_queries = get_parquet_row_counts(partial_path)
                partial_cache_path = partial_path
                log_info(f"Partial cache found: {p_refs} refs, {p_queries} queries. "
                         f"Will resume loading remaining spectra.", cfg.quiet)

    # =========================================================================
    # Phase 2: Peptide search (if not Parquet)
    # =========================================================================
    search_results = None

    if not use_parquet:
        # Check for search checkpoint first
        search_checkpoint_loaded = False
        if cfg.cache_enabled and cache_fingerprint:
            from src.parquet_cache import load_search_checkpoint
            checkpoint_results = load_search_checkpoint(cache_dir, cache_fingerprint)
            if checkpoint_results:
                search_results = checkpoint_results
                search_checkpoint_loaded = True
                total_psms = sum(len(vg.matches) for vg in search_results.values())
                log_info(f"Search checkpoint loaded: {total_psms} PSMs across "
                         f"{len(search_results)} variants (skipping idXML search)", cfg.quiet)

        if not search_checkpoint_loaded:
            from src.peptide_search import search_peptides_in_idxmls

            log_info("Searching peptides in idXML files...", cfg.quiet)
            search_errors = []
            search_progress_bar = TerminalProgress(
                len(cfg.idxml_files), "Search", quiet=cfg.quiet
            )

            def search_progress(current, total, message):
                search_progress_bar.update(current)

            def search_error(msg):
                search_errors.append(msg)

            search_results = search_peptides_in_idxmls(
                cfg.idxml_files,
                cfg.peptides,
                mzml_paths=cfg.mzml_files,
                progress_callback=search_progress,
                error_callback=search_error,
                workers=cfg.workers,
            )
            search_progress_bar.finish("Search complete")

            if not search_results:
                log_error("No peptide matches found in idXML files")
                return 1

            total_psms = sum(len(vg.matches) for vg in search_results.values())
            log_info(f"Found {total_psms} PSMs across {len(search_results)} peptide variants", cfg.quiet)

            if search_errors and cfg.verbose:
                for err in search_errors[:5]:
                    log_debug(err, True)
                if len(search_errors) > 5:
                    log_debug(f"... and {len(search_errors) - 5} more warnings", True)

            # Save search checkpoint immediately
            if cfg.cache_enabled and cache_fingerprint:
                from src.parquet_cache import save_search_checkpoint
                saved = save_search_checkpoint(cache_dir, cache_fingerprint, search_results)
                if saved:
                    log_info(f"Search checkpoint saved: {Path(saved).name}", cfg.quiet)

    # =========================================================================
    # Phase 3: Prepare for processing
    # =========================================================================
    if not use_parquet:
        variant_groups = copy.deepcopy(search_results)
        mzml_paths = [Path(m) for m in cfg.mzml_files]
        total_psms = sum(len(vg.matches) for vg in variant_groups.values())
    else:
        variant_groups = None
        mzml_paths = None
        total_psms = 0

    # Build run settings dict (everything that affects output)
    from src.export import find_or_create_export_directory

    run_settings = {
        "idxml_files": sorted(cfg.idxml_files),
        "mzml_files": sorted(cfg.mzml_files),
        "peptides": sorted(p.upper() for p in cfg.peptides),
        "reference_substring": cfg.reference_substring,
        "ion_types": sorted(cfg.ion_types),
        "mz_tolerance": cfg.mz_tolerance,
        "add_losses": cfg.add_losses,
        "add_first_prefix_ion": cfg.add_first_prefix_ion,
        "remove_precursor_peaks": cfg.remove_precursor_peaks,
        "intensity_mode": cfg.intensity_mode,
        "export_formats": sorted(cfg.export_formats),
        "dpi": cfg.dpi,
        "figsize": list(cfg.figsize),
        "plot_unpaired": cfg.plot_unpaired,
    }

    # Numbered folder: ProCa_2_run_001/, _002/, etc. Reuses matching folder.
    root_dir, png_dir, svg_dir = find_or_create_export_directory(
        cfg.export_dir, cfg.project_name, run_settings
    )
    log_info(f"Export directory: {root_dir}", cfg.quiet)

    # Count expected refs/queries for progress
    if use_parquet:
        from src.parquet_cache import get_parquet_row_counts
        expected_refs, expected_queries = get_parquet_row_counts(cfg.parquet_file)
    else:
        ref_names_lower = {Path(m).name.lower() for m in cfg.mzml_files
                           if cfg.reference_substring and cfg.reference_substring.lower() in Path(m).name.lower()}
        ref_paths = {str(Path(m).resolve()).lower() for m in cfg.mzml_files
                     if cfg.reference_substring and cfg.reference_substring.lower() in Path(m).name.lower()}
        expected_refs = sum(
            1 for vg in variant_groups.values() for m in vg.matches
            if (m.mzml_resolved_path and str(Path(m.mzml_resolved_path).resolve()).lower() in ref_paths)
            or (not m.mzml_resolved_path and m.mzml_filename and m.mzml_filename.lower() in ref_names_lower)
        )
        expected_queries = total_psms - expected_refs

    if expected_refs > 0:
        log_info(f"Expected: {expected_refs} reference + {expected_queries} query spectra", cfg.quiet)
    else:
        log_info(f"Expected: {expected_queries} spectra (no reference)", cfg.quiet)

    # =========================================================================
    # Phase 4: Process spectra (streaming)
    # =========================================================================
    start_time = time.time()
    ref_start_time = start_time
    query_start_time = [None]  # mutable for closure

    # Progress trackers
    ref_progress = TerminalProgress(expected_refs, "References", quiet=cfg.quiet) if expected_refs > 0 else None
    query_progress = TerminalProgress(
        expected_queries,
        "Queries" if expected_refs > 0 else "Spectra",
        quiet=cfg.quiet,
    )

    # Results tracking
    results_summary = {
        "reference_spectra": 0,
        "query_spectra": 0,
        "plots_generated": 0,
        "files_exported": {"png": [], "svg": []},
        "errors": [],
        "export_root": str(root_dir),
    }

    annotated_references: Dict[str, list] = {}
    refs_done = [expected_refs == 0]  # mutable for closure
    summary_rows: List[dict] = []
    _cached_queries: List[Tuple[str, object]] = []
    _cached_references: Dict[str, list] = {}  # accumulated during streaming for cache saves

    # Task collection for parallel processing
    from src.parallel_processor import ComparisonTask
    comparison_tasks: List[ComparisonTask] = []
    task_counter = [0]

    # --- Annotation helper ---
    def annotate_spectrum(spectrum, is_reference=False):
        spec_type = "reference" if is_reference else "query"

        if spectrum.mz is None or len(spectrum.mz) == 0:
            results_summary["errors"].append(f"Empty m/z: {spec_type} {spectrum.peptidoform}")
            return None
        if spectrum.intensity is None or len(spectrum.intensity) == 0:
            results_summary["errors"].append(f"Empty intensity: {spec_type} {spectrum.peptidoform}")
            return None
        if not spectrum.peptidoform:
            results_summary["errors"].append(f"Missing peptidoform: {spec_type}")
            return None
        if spectrum.precursor_charge <= 0:
            results_summary["errors"].append(f"Invalid charge: {spec_type} {spectrum.peptidoform}")
            return None

        try:
            return prepare_spectrum_for_plot(
                mz=spectrum.mz,
                intensity=spectrum.intensity,
                peptidoform=spectrum.peptidoform,
                charge=spectrum.precursor_charge,
                ion_types=cfg.ion_types,
                mz_tolerance=cfg.mz_tolerance,
                add_losses=cfg.add_losses,
                add_first_prefix_ion=cfg.add_first_prefix_ion,
                remove_precursor_peaks=cfg.remove_precursor_peaks,
                normalize=cfg.normalize,
            )
        except Exception as e:
            results_summary["errors"].append(f"Annotation error ({spec_type} {spectrum.peptidoform}): {e}")
            return None

    # --- Cache save helper (supports partial + complete saves) ---
    _cache_complete = [False]  # True once a complete (non-partial) save succeeds

    def _try_save_cache(stage: str, partial: bool = False):
        """Save Parquet cache. Overwrites previous saves with updated data.

        Args:
            stage: Label for logging (e.g. "after refs", "incremental 500")
            partial: If True, marks cache as incomplete (more queries expected)
        """
        if _cache_complete[0]:
            return  # Complete cache already saved, no need to overwrite
        if (loaded_from_cache
                or use_parquet
                or not cfg.cache_enabled
                or not cache_fingerprint):
            return
        # For partial saves, need at least references; for complete, need queries
        has_refs = _cached_references and any(_cached_references.values())
        if not has_refs and not _cached_queries:
            return
        try:
            from src.parquet_cache import save_spectra_to_parquet

            cache_path = str(Path(cache_dir) / f"{cache_fingerprint}.parquet")
            cache_metadata = {
                "fingerprint": cache_fingerprint,
                "export_root": str(root_dir),
                "idxml_files": cfg.idxml_files,
                "mzml_files": cfg.mzml_files,
                "target_peptides": cfg.peptides,
                "reference_substring": cfg.reference_substring,
            }
            saved = save_spectra_to_parquet(
                reference_spectra=_cached_references,
                query_spectra=_cached_queries,
                output_path=cache_path,
                cache_metadata=cache_metadata,
                partial=partial,
            )
            if saved:
                if not partial:
                    _cache_complete[0] = True
                n_refs = sum(len(v) for v in _cached_references.values())
                tag = "partial" if partial else "complete"
                log_info(f"Cache saved [{tag}] ({stage}): {Path(saved).name} "
                         f"({len(_cached_queries)} queries, {n_refs} refs)",
                         cfg.quiet)
        except Exception as cache_err:
            log_info(f"Warning: cache save failed ({stage}): {cache_err}", cfg.quiet)

    # --- Callbacks ---
    # Incremental cache save interval
    _INCREMENTAL_SAVE_INTERVAL = 500
    _last_incremental_save = [0]  # query count at last incremental save

    def on_reference_loaded(variant_key: str, count: int):
        results_summary["reference_spectra"] += count
        if ref_progress:
            ref_progress.update(results_summary["reference_spectra"])

        if results_summary["reference_spectra"] >= expected_refs and not refs_done[0]:
            refs_done[0] = True
            if ref_progress:
                ref_progress.finish(f"References: {results_summary['reference_spectra']} loaded")
            query_start_time[0] = time.time()

    def on_query_ready(result: StreamingResult):
        """Annotate query + build ComparisonTasks (no figure creation here)."""
        query = result.query_spectrum
        refs = result.reference_spectra

        # Accumulate references into cache (once per variant)
        if not use_parquet and result.variant_key not in _cached_references:
            _cached_references[result.variant_key] = list(refs)

        # Accumulate queries for cache + incremental save every N queries
        if not use_parquet:
            _cached_queries.append((result.variant_key, query))
            n_queries = len(_cached_queries)
            if (n_queries - _last_incremental_save[0]) >= _INCREMENTAL_SAVE_INTERVAL:
                _last_incremental_save[0] = n_queries
                _try_save_cache(f"incremental {n_queries}", partial=True)

        # Mark refs done on first query
        if not refs_done[0]:
            refs_done[0] = True
            if ref_progress:
                ref_progress.finish(f"References: {results_summary['reference_spectra']} loaded")
            query_start_time[0] = time.time()

        results_summary["query_spectra"] += 1
        query_progress.update()

        # Annotate query
        query_df = annotate_spectrum(query, is_reference=False)
        if query_df is None:
            return

        # Ensure references are annotated (cached per variant)
        if result.variant_key not in annotated_references:
            annotated_references[result.variant_key] = []
            for ref in refs:
                ref_df = annotate_spectrum(ref, is_reference=True)
                if ref_df is not None:
                    annotated_references[result.variant_key].append((ref, ref_df))

        ref_data = annotated_references.get(result.variant_key, [])

        # Folder and title
        if result.modification_key:
            folder_name = f"{result.base_sequence}_{sanitize_filename(result.modification_key)}"
            plot_title = f"{result.base_sequence} [{result.modification_key}]"
        else:
            folder_name = result.base_sequence
            plot_title = result.base_sequence

        if not ref_data:
            if not cfg.plot_unpaired:
                log_debug(f"  Skipped (no reference): {Path(query.source_file).stem}", cfg.verbose)
                return

            # Build single spectrum task
            query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
            single_label = (f"{Path(query.source_file).stem}"
                            f"  (RT: {query_rt}, ID: {query.spectrum_id})")
            query_stem = Path(query.source_file).stem[:30]
            base_filename = f"{result.base_sequence}_{query_stem}"
            safe_filename = sanitize_filename(base_filename)

            comparison_tasks.append(ComparisonTask(
                task_id=task_counter[0],
                task_type="single",
                query_df=query_df,
                query_label=single_label,
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=str(png_dir),
                svg_dir=str(svg_dir),
                export_formats=cfg.export_formats,
                figsize=cfg.figsize,
                dpi=cfg.dpi,
                normalize=cfg.normalize,
                query_file=Path(query.source_file).name,
                peptide=result.base_sequence,
                modification=result.modification_key or "",
                charge=query.precursor_charge,
            ))
            task_counter[0] += 1
            return

        # Build mirror task for each reference
        for ref, ref_df in ref_data:
            query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
            ref_rt = f"{ref.rt:.1f}" if ref.rt else "N/A"
            query_label = (f"Query: {Path(query.source_file).stem}"
                           f"  (RT: {query_rt}, ID: {query.spectrum_id})")
            ref_label = (f"Ref: {Path(ref.source_file).stem}"
                         f"  (RT: {ref_rt}, ID: {ref.spectrum_id})")

            query_stem = Path(query.source_file).stem[:30]
            ref_stem = Path(ref.source_file).stem[:30]
            base_filename = f"{result.base_sequence}_{query_stem}_vs_{ref_stem}"
            safe_filename = sanitize_filename(base_filename)

            comparison_tasks.append(ComparisonTask(
                task_id=task_counter[0],
                task_type="mirror",
                query_df=query_df,
                ref_df=ref_df,
                query_label=query_label,
                ref_label=ref_label,
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=str(png_dir),
                svg_dir=str(svg_dir),
                export_formats=cfg.export_formats,
                figsize=cfg.figsize,
                dpi=cfg.dpi,
                normalize=cfg.normalize,
                query_file=Path(query.source_file).name,
                ref_file=Path(ref.source_file).name,
                peptide=result.base_sequence,
                modification=result.modification_key or "",
                charge=query.precursor_charge,
            ))
            task_counter[0] += 1

    def on_file_start(filename: str, is_reference: bool):
        file_type = "reference" if is_reference else "query"
        log_debug(f"Opening {file_type} file: {filename}", cfg.verbose)

    def on_file_complete(filename: str, loaded: int, not_found: int):
        log_debug(f"Completed {filename}: {loaded} loaded, {not_found} not found", cfg.verbose)

    def on_error(error_msg: str):
        results_summary["errors"].append(error_msg)
        log_debug(f"Error: {error_msg}", cfg.verbose)

    # --- Phase 4a: Stream spectra and annotate (collect tasks) ---
    log_info("", cfg.quiet)  # blank line before progress

    try:
        if use_parquet:
            log_info(f"Loading from Parquet: {Path(cfg.parquet_file).name}", cfg.quiet)
            _parquet_refs = load_spectra_from_parquet(
                parquet_path=cfg.parquet_file,
                reference_substring=cfg.reference_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                error_callback=on_error,
            )
            _cached_references.update(_parquet_refs)
        elif partial_cache_path:
            # Resume from partial cache: load cached data, then stream remaining
            log_info(f"Resuming from partial cache: {Path(partial_cache_path).name}", cfg.quiet)

            # Step 1: Load partial Parquet through callbacks (populates tasks/annotations)
            _parquet_refs = load_spectra_from_parquet(
                parquet_path=partial_cache_path,
                reference_substring=cfg.reference_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                error_callback=on_error,
            )
            _cached_references.update(_parquet_refs)

            cached_refs = results_summary["reference_spectra"]
            cached_queries = results_summary["query_spectra"]
            log_info(f"Partial cache loaded: {cached_refs} refs, {cached_queries} queries", cfg.quiet)

            # Step 2: Build skip set from partial cache
            from src.parquet_cache import get_cached_spectrum_ids
            skip_ids = get_cached_spectrum_ids(partial_cache_path)
            log_info(f"Skip set: {len(skip_ids)} spectra already cached", cfg.quiet)

            # Step 3: Stream remaining spectra from mzML
            _streaming_refs = load_spectra_streaming(
                variant_groups=variant_groups,
                mzml_paths=mzml_paths,
                reference_substring=cfg.reference_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                on_file_start=on_file_start,
                on_file_complete=on_file_complete,
                error_callback=on_error,
                workers=cfg.workers,
                skip_ids=skip_ids,
            )
            # Merge any new references
            for vk, specs in _streaming_refs.items():
                if specs:
                    if vk not in _cached_references:
                        _cached_references[vk] = []
                    _cached_references[vk].extend(specs)
        else:
            _streaming_refs = load_spectra_streaming(
                variant_groups=variant_groups,
                mzml_paths=mzml_paths,
                reference_substring=cfg.reference_substring,
                on_reference_loaded=on_reference_loaded,
                on_query_ready=on_query_ready,
                on_file_start=on_file_start,
                on_file_complete=on_file_complete,
                error_callback=on_error,
                workers=cfg.workers,
            )

        # Finish progress bars
        if ref_progress and not refs_done[0]:
            ref_progress.finish(f"References: {results_summary['reference_spectra']} loaded")
        query_progress.finish(
            f"{'Queries' if expected_refs > 0 else 'Spectra'}: "
            f"{results_summary['query_spectra']} annotated, "
            f"{len(comparison_tasks)} plots queued"
        )

    except Exception as e:
        log_error(f"Processing failed: {e}")
        logger.exception("Processing error")
        return 1

    # Save point: after streaming complete — final complete cache
    _try_save_cache("after loading", partial=False)

    # --- Phase 4b: Generate plots in parallel ---
    if comparison_tasks:
        from src.parallel_processor import run_parallel_comparisons
        from src.tmpfs_staging import TmpfsStaging, StagingConfig

        # Check for existing plot files (resumable runs)
        # Uses final disk paths — before any tmpfs rewriting
        skipped = 0
        for task in comparison_tasks:
            outputs_exist = True
            if "png" in task.export_formats and task.png_dir:
                if not (Path(task.png_dir) / task.folder_name / f"{task.safe_filename}.png").exists():
                    outputs_exist = False
            if "svg" in task.export_formats and task.svg_dir:
                if not (Path(task.svg_dir) / task.folder_name / f"{task.safe_filename}.svg").exists():
                    outputs_exist = False
            if outputs_exist:
                task.skip_plot = True
                skipped += 1

        if skipped > 0:
            log_info(f"Skipping {skipped}/{len(comparison_tasks)} plots (already on disk)", cfg.quiet)

        pending = len(comparison_tasks) - skipped
        log_info("", cfg.quiet)
        log_info(f"Generating {pending} plots, {skipped} cached ({cfg.workers} workers)...", cfg.quiet)

        plot_progress = TerminalProgress(len(comparison_tasks), "Plots", quiet=cfg.quiet)

        # Global counter for progress across chunks
        _global_completed = [0]

        def on_plot_complete_chunked(result, completed, total):
            _global_completed[0] += 1
            plot_progress.update(_global_completed[0])
            if result.success:
                results_summary["plots_generated"] += 1
                if result.png_path:
                    results_summary["files_exported"]["png"].append(result.png_path)
                if result.svg_path:
                    results_summary["files_exported"]["svg"].append(result.svg_path)
                if result.summary_row:
                    summary_rows.append(result.summary_row)
            else:
                results_summary["errors"].append(f"Plot error: {result.error}")

        staging_config = StagingConfig(
            enabled=cfg.tmpfs_staging,
            tmpfs_base=cfg.tmpfs_path,
            ram_fraction=cfg.tmpfs_ram_fraction,
        )

        with TmpfsStaging(str(png_dir), str(svg_dir), cfg.export_formats,
                          staging_config,
                          log_fn=lambda msg: log_info(msg, cfg.quiet)) as staging:
            for chunk in staging.iter_chunks(comparison_tasks):
                run_parallel_comparisons(
                    tasks=chunk,
                    workers=cfg.workers,
                    on_task_complete=on_plot_complete_chunked,
                )
                staging.flush()

        plot_progress.finish(
            f"Plots: {results_summary['plots_generated']} generated"
        )

    # Save point 2: after plots
    _try_save_cache("after plots")

    total_time = time.time() - start_time

    # =========================================================================
    # Phase 5: Write comparison summary CSV + XLSX
    # =========================================================================
    if summary_rows:
        try:
            import pandas as _pd
            summary_df = _pd.DataFrame(summary_rows)

            all_filenames = list(set(
                summary_df["query_file"].tolist() + summary_df["ref_file"].tolist()
            ))
            short_map = dict(zip(all_filenames, shorten_filenames(all_filenames)))
            summary_df.insert(0, "query_short", summary_df["query_file"].map(short_map))
            summary_df.insert(1, "ref_short", summary_df["ref_file"].map(short_map))

            col_order = [
                "query_short", "ref_short", "peptide", "modification", "charge",
                "R2m", "R2b", "Cos", "Ent", "Dice", "Jacc", "ShR",
                "total_q", "total_r", "ann_q", "ann_r", "shared",
                "query_file", "ref_file",
            ]
            summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]

            rename_map = {
                "query_short": "Query", "ref_short": "Reference",
                "peptide": "Peptide", "modification": "Modification",
                "charge": "Charge",
                "R2m": "R² match", "R2b": "R² binned",
                "Cos": "Cosine", "Ent": "Entropy",
                "Dice": "Dice", "Jacc": "Jaccard", "ShR": "SharedRatio",
                "total_q": "Peaks (Q)", "total_r": "Peaks (R)",
                "ann_q": "Annotated (Q)", "ann_r": "Annotated (R)",
                "shared": "Shared Ions",
                "query_file": "Query File", "ref_file": "Ref File",
            }
            summary_df = summary_df.rename(columns=rename_map)

            csv_path = root_dir / "comparison_summary.csv"
            summary_df.to_csv(csv_path, index=False)
            results_summary["summary_csv"] = str(csv_path)
            log_info(f"Summary CSV: {csv_path}", cfg.quiet)

            try:
                xlsx_path = root_dir / "comparison_summary.xlsx"
                summary_df.to_excel(xlsx_path, index=False, sheet_name="Comparisons")
                results_summary["summary_xlsx"] = str(xlsx_path)
                log_info(f"Summary XLSX: {xlsx_path}", cfg.quiet)
            except Exception as xlsx_err:
                log_debug(f"Could not write XLSX: {xlsx_err}", cfg.verbose)

        except Exception as summary_err:
            log_debug(f"Could not write summary: {summary_err}", cfg.verbose)

    # Save point 3: after summary export
    _try_save_cache("after summary")

    # =========================================================================
    # Final summary
    # =========================================================================
    log_info("", cfg.quiet)
    log_info("=" * 60, cfg.quiet)
    log_info("COMPLETE", cfg.quiet)
    log_info("=" * 60, cfg.quiet)
    log_info(f"References:  {results_summary['reference_spectra']}", cfg.quiet)
    log_info(f"Queries:     {results_summary['query_spectra']}", cfg.quiet)
    log_info(f"Plots:       {results_summary['plots_generated']}", cfg.quiet)
    log_info(f"PNG files:   {len(results_summary['files_exported']['png'])}", cfg.quiet)
    log_info(f"SVG files:   {len(results_summary['files_exported']['svg'])}", cfg.quiet)
    log_info(f"Total time:  {format_time(total_time)}", cfg.quiet)
    log_info(f"Output:      {root_dir}", cfg.quiet)

    if results_summary["errors"]:
        n_err = len(results_summary["errors"])
        log_info(f"Warnings:    {n_err}", cfg.quiet)
        if cfg.verbose:
            for err in results_summary["errors"][:10]:
                log_debug(err, True)
            if n_err > 10:
                log_debug(f"... and {n_err - 10} more", True)

    if results_summary["plots_generated"] == 0:
        log_error("No plots were generated. Check your inputs and reference substring.")
        return 1

    return 0


# ============================================================================
# Argument parser + main
# ============================================================================

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mspecviz-cli",
        description="MSpecViz CLI — Terminal-based MS/MS spectrum analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python cli.py --config config.yaml
  python cli.py --config config.yaml --peptides peptides.txt --verbose
  python cli.py --config config.yaml --mz-tolerance 0.05 --add-losses --quiet
""",
    )

    # Required
    parser.add_argument("--config", required=True, help="Path to YAML configuration file")

    # Input overrides
    parser.add_argument("--idxml", nargs="+", help="idXML file paths (overrides config)")
    parser.add_argument("--mzml", nargs="+", help="mzML file paths (overrides config)")
    parser.add_argument("--parquet", help="Parquet file path (skips search + mzML loading)")
    parser.add_argument("--peptides", help="Path to peptides .txt file (one per line, overrides config)")
    parser.add_argument("--reference-substring", help="Reference file identifier substring")

    # Annotation overrides
    parser.add_argument("--ion-types", nargs="+", help="Ion types to annotate (e.g., b y)")
    parser.add_argument("--mz-tolerance", type=float, help="m/z tolerance in Da")
    parser.add_argument("--add-losses", action="store_true", help="Add neutral losses (H2O, NH3)")
    parser.add_argument("--no-first-prefix-ion", action="store_true", help="Disable first prefix ion")
    parser.add_argument("--remove-precursor-peaks", action="store_true", help="Remove precursor peaks")
    parser.add_argument("--raw-intensity", action="store_true", help="Use raw intensities (not normalized)")

    # Export overrides
    parser.add_argument("--export-dir", help="Export directory path")
    parser.add_argument("--project-name", help="Project name for export folder")
    parser.add_argument("--formats", nargs="+", choices=["png", "svg"], help="Export formats")
    parser.add_argument("--dpi", type=int, help="DPI for image export")
    parser.add_argument("--plot-unpaired", action="store_true", help="Plot spectra without matching reference")

    # Cache
    parser.add_argument("--no-cache", action="store_true", help="Disable Parquet auto-caching")

    # Performance
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel worker processes for plot generation (default: 4)")
    parser.add_argument("--no-tmpfs", action="store_true",
                        help="Disable RAM-backed tmpfs staging for plot writes")

    # Verbosity
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose/debug output")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress all output except errors")

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    # Validate config file exists
    if not Path(args.config).exists():
        log_error(f"Config file not found: {args.config}")
        sys.exit(2)

    # Setup logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s %(name)s: %(message)s")
    elif args.quiet:
        logging.basicConfig(level=logging.ERROR, format="%(levelname)s: %(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    # Load and merge configuration
    try:
        cfg = load_config(args.config, args)
    except Exception as e:
        log_error(f"Configuration error: {e}")
        sys.exit(2)

    # Run pipeline
    exit_code = run_pipeline(cfg)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
