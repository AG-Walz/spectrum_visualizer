"""Shared export orchestration (called by cli.main and streamlit_app.views.export).
- `run_settings.json` schematizes and 
- Export-folder allocation/renaming with settings fingerprinting
  (`find_or_create_export_directory`, `commit_export_directory`,
  `sanitize_filename`); folders follow {project}_{NNN}_{YYYYMMDD}/ and reuse
  existing folders when the settings fingerprint matches.
- The top-level export pipeline (`run_export` + the `ExportSelection`,
  `PlotSettings`, `ExportCallbacks`, `ExportResult` dataclasses) that writes
  parquet/CSV/XLSX, dispatches mirror-plot rendering to a worker pool, and
  emits the summary plots.

Dispatches to core.pipeline.build_plot_tasks and
core.plotting.process_single_task / generate_all_summary_plots.
"""
from __future__ import annotations

import hashlib
import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from core.defaults import DEFAULTS
from core.helpers import run_tasks
from core.pipeline import build_plot_tasks
from core.plotting import generate_all_summary_plots, process_single_task

logger = logging.getLogger(__name__)


# ============================================================================
# run_settings.json schema + builder
# ============================================================================

# A frontend's parameter source, adapted to a uniform (key, default) -> value
# lookup — the CLI passes getattr(config, ...), Streamlit reads params.
ParamGetter = Callable[[str, Any], Any]


# Each entry: (run_settings_key, default). Frontends map their internal storage
# (CLI SimpleNamespace attr, Streamlit params widget key) onto these keys
# via their getter — see cli.main and streamlit_app.utils.page_state.
RUN_SETTINGS_FIELDS: list[tuple[str, object]] = [
    ("idxml_files",            []),
    ("mzml_files",             []),
    ("reference_substring",    ""),
    ("ion_types",              list(DEFAULTS.ion_types)),
    ("mz_tolerance",           DEFAULTS.mz_tolerance),
    ("add_neutral_losses",     DEFAULTS.add_neutral_losses),
    ("intensity_mode",         DEFAULTS.intensity_mode),
    ("koina_enabled",          DEFAULTS.koina_enabled),
    ("koina_model",            DEFAULTS.koina_model),
    ("koina_collision_energy", DEFAULTS.koina_collision_energy),
    ("koina_instrument_type",  DEFAULTS.koina_instrument_type),
    ("export_formats",         list(DEFAULTS.export_formats)),
    ("export_resolution",      DEFAULTS.export_resolution),
    ("figure_size",            list(DEFAULTS.figure_size)),
    ("summary_metric",         DEFAULTS.summary_metric),
    ("summary_plot_types",     list(DEFAULTS.summary_plot_types)),
]


def build_run_settings(getter: ParamGetter, peptides: list[str]) -> dict:
    """Assemble the run_settings.json dict from any frontend's parameter source.
    peptides is passed separately because frontends parse it differently (YAML list vs newline string).
    The returned dict holds every RUN_SETTINGS_FIELDS key plus peptides and the
    underscore-prefixed timestamps _created_at / _last_run (both set to now on creation).
    """
    # Pull each declared field through the frontend getter, falling back to its default.
    settings: dict = {key: getter(key, default) for key, default in RUN_SETTINGS_FIELDS}
    settings["peptides"] = peptides
    # Both timestamps start equal; _last_run is later refreshed when an export reuses this folder.
    now = datetime.now().isoformat()
    settings["_created_at"] = now
    settings["_last_run"] = now
    return settings


# ============================================================================
# Export-folder allocation (fingerprinted versioned dirs)
# ============================================================================

def sanitize_filename(name: str) -> str:
    """Make user provided name safe to use as a filesystem path component.
    """
    # Replace each filesystem-unsafe character with an underscore
    sanitized = re.sub(r'[<>:"/\\|?*\[\]]', '_', name)
    # Collapse consecutive underscores, then strip stray edge underscores/spaces
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_ ')
    # Return cleaned string with a maximum length of 100 characters
    return sanitized[:100]


def _settings_fingerprint(run_settings: dict) -> str:
    """Return a 16-character SHA-256 hash identifying a settings dictionary.
    Identical settings always yield the same hash regardless of key order.
    Used to detect re-runs with unchanged settings.
    """
    # Canonical JSON: sorted keys make the hash order-independent
    canonical = json.dumps(run_settings, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def _entry_fingerprint(settings_file: Path) -> str | None:
    """Read the _fingerprint value stored in a run_settings.json file.
    Returns the fingerprint string, or None if the file is missing.
    """
    try:
        with open(settings_file, encoding="utf-8") as f:
            return json.load(f).get("_fingerprint")
    except (json.JSONDecodeError, OSError) as error:
        logger.warning(f"Could not read {settings_file}: {error}")
        return None


def _scan_existing_dirs(
    base: Path, pattern: re.Pattern, fingerprint: str,
) -> tuple[int, Path | None, int | None]:
    """Scan base for existing export folders matching the project naming scheme.
    Walks the immediate children of the user-defined base folder, considering only directories whose
    name match the regex pattern {project}_{NNN}_{YYYYMMDD}. Returns max_number (highest version number seen),
    match_entry and match_number point matching the stored fingerprint.
    """
    max_number = 0
    match_entry = None
    match_number = None
    # A missing base directory means no prior runs — nothing to scan
    if not base.exists():
        return max_number, match_entry, match_number

    # Iterate sorted children so the fingerprint match is deterministic
    for entry in sorted(base.iterdir()):
        if not entry.is_dir():
            continue
        # Skip directories that do not follow the versioned naming scheme
        name_match = pattern.match(entry.name)
        if not name_match:
            continue
        # Track the highest version number for the fallback allocation
        number = int(name_match.group(1))
        max_number = max(max_number, number)
        # Record the first folder whose settings fingerprint matches this run
        if match_entry is None:
            settings_file = entry / "run_settings.json"
            if settings_file.exists() and _entry_fingerprint(settings_file) == fingerprint:
                match_entry = entry
                match_number = number
    return max_number, match_entry, match_number


def _rename_or_reuse_match(prefix: str, match_entry: Path,
                           match_number: int, today: str) -> Path:
    """Re-date a fingerprint-matching export folder to current (today's) date while keeping its version number.
    """
    new_path = match_entry.with_name(f"{prefix}_{match_number:03d}_{today}")
    # Folder is already dated today — nothing to rename
    if match_entry == new_path:
        logger.info(f"Reusing export directory: {new_path.name}")
        return new_path
    # Re-date the folder; on failure keep using the original (rename is not critical)
    try:
        match_entry.rename(new_path)
        logger.info(f"Renamed {match_entry.name} -> {new_path.name}")
        return new_path
    except OSError as error:
        logger.warning(f"Could not rename folder: {error}; reusing original: {match_entry.name}")
        return match_entry


def find_or_create_export_directory(
    base_dir: str,
    project_name: str,
    run_settings: dict,
) -> Path:
    """Decide the export folder for a run.
    Returns a Path named {project}_{NNN}_{YYYYMMDD}. When a prior folder
    carries the same settings fingerprint, that folder is reused (re-dated to today
    via _rename_or_reuse_match); otherwise the next free version number
    max_number + 1 is allocated. This function only plans the path — it does not
    create the folder on disk; call commit_export_directory for that.
    """
    # Setup: sanitize project name, compute fingerprint, build name regex
    base = Path(base_dir)
    prefix = sanitize_filename(project_name)
    fingerprint = _settings_fingerprint(run_settings)
    today = datetime.now().strftime("%Y%m%d")
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d{{3}})_(\d{{8}})$")

    # Scan: highest existing NNN and any fingerprint-matching folder
    max_number, match_entry, match_number = _scan_existing_dirs(base, pattern, fingerprint)

    # Reuse on match, otherwise allocate next NNN
    if match_entry is not None:
        return _rename_or_reuse_match(prefix, match_entry, match_number, today)

    next_number = max_number + 1
    return base / f"{prefix}_{next_number:03d}_{today}"


def commit_export_directory(root_dir: Path, run_settings: dict) -> None:
    """Materialise the export folder on disk and persist its run_settings.json.
    Creates root_dir (and parents) if needed. For a freshly allocated folder it
    writes a new run_settings.json stamped with the settings fingerprint and a
    creation timestamp; for a reused folder it leaves the existing settings intact
    and only refreshes the _last_run timestamp. Returns nothing — I/O failures
    are logged as warnings rather than raised, so a failed write never aborts the export.
    """
    # Ensure the export folder exists
    root_dir.mkdir(parents=True, exist_ok=True)

    # Reused folder: bump _last_run on the existing settings file
    settings_file = root_dir / "run_settings.json"
    if settings_file.exists():
        try:
            with open(settings_file, encoding="utf-8") as f:
                existing_settings = json.load(f)
            existing_settings["_last_run"] = datetime.now().isoformat()
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(existing_settings, f, indent=2, default=str)
        except (json.JSONDecodeError, OSError) as error:
            logger.warning(f"Could not update run_settings.json: {error}")
    # Fresh folder: stamp fingerprint + creation time, write new file
    else:
        settings_to_save = {
            "_fingerprint": _settings_fingerprint(run_settings),
            "_created_at": datetime.now().isoformat(),
            **run_settings,
        }
        try:
            with open(settings_file, "w", encoding="utf-8") as f:
                json.dump(settings_to_save, f, indent=2, default=str)
        except OSError as error:
            logger.warning(f"Could not write run_settings.json: {error}")


# ============================================================================
# Export pipeline (run_export + request/result dataclasses)
# ============================================================================

# Two shapes a frontend may hand in: a list [(psm_id, spec, df), ...] or a
# dict {psm_id: (spec, df)} — run_export normalises both into a list.
AnnotatedQueries = list[tuple] | dict[str, tuple]
# The four openms_insight LazyFrames written one parquet file each.
_PARQUET_FRAMES = ("psm", "peaks", "sequence", "heatmap")
# Raster/vector image formats accepted for plots — distinct from html,
# which is handled on its own track because it is interactive, not a picture.
_PICTURE_FORMATS = ("png", "svg", "pdf")


@dataclass
class ExportSelection:
    """Which artefact groups to write (parquet / CSV / mirror plots / summaries / HTML).
    """
    parquet: bool = True
    csv: bool = True
    mirror_plots: bool = True
    summary_plots: bool = True
    html: bool = False


@dataclass
class PlotSettings:
    """Figure size, image formats, DPI, and plot-worker count.
    """
    export_formats: list[str] = field(default_factory=lambda: ["png", "svg"])
    figure_size: tuple[int, int] = (12, 8)
    export_resolution: int = 300
    plot_unpaired_queries: bool = False
    workers: int = 1


@dataclass
class ExportCallbacks:
    """Optional UI hooks for step messages, plot progress, and warnings (no-op when None).
    """
    on_step: Callable[[str], None] | None = None
    on_plot_progress: Callable[[int, int], None] | None = None
    on_warning: Callable[[str], None] | None = None


@dataclass
class ExportResult:
    """Final outcome: the chosen folder, every saved file, and per-group plot counts.
    """
    export_root: Path
    saved_files: list[str]
    mirror_plot_count: int = 0
    summary_plot_count: int = 0


def _resolve_format_dirs(
    root_dir: Path, export_formats: list[str],
) -> tuple[dict[str, str], str]:
    """Map each requested format to its subfolder path under root_dir.
    Returns (format_dirs, html_dir) where format_dirs is {fmt: path} for
    every picture format present in export_formats and html_dir is the html
    subfolder path (empty string when html is not requested). Paths are strings only —
    the directories themselves are created later by the plot workers.
    """
    format_dirs = {fmt: str(root_dir / fmt) for fmt in _PICTURE_FORMATS if fmt in export_formats}
    html_dir = str(root_dir / "html") if "html" in export_formats else ""
    return format_dirs, html_dir


def _run_plot_tasks(
    annotated_queries: list[tuple],
    reference_store: dict,
    variant_groups: dict,
    root_dir: Path,
    plot_settings: PlotSettings,
    do_html: bool,
    saved_files: list[str],
    on_progress: Callable[[int, int], None],
    on_warning: Callable[[str], None],
) -> int:
    """Build one mirror plot per (query, reference) pair and render via a worker pool.
    Each successfully rendered task contributes its non-empty per-format output paths
    to saved_files (mutated in place). Returns the count of files written; 0 when
    no tasks were produced. Worker failures are reported through on_warning rather
    than raised, so a single bad PSM does not abort the whole export.
    """
    # Resolve format set: append html when requested, map each to a subdir
    base_formats = list(plot_settings.export_formats)
    if do_html and "html" not in base_formats:
        base_formats.append("html")
    format_dirs, html_dir = _resolve_format_dirs(root_dir, base_formats)

    # Build tasks: one per (query, reference) pair with all per-format output paths
    tasks = build_plot_tasks(
        annotated_queries, reference_store, variant_groups,
        sanitize_fn=sanitize_filename,
        format_dirs=format_dirs,
        figure_size=plot_settings.figure_size,
        export_resolution=plot_settings.export_resolution,
        html_dir=html_dir,
        plot_unpaired_queries=plot_settings.plot_unpaired_queries,
    )
    if not tasks:
        return 0

    # Dispatch: tally saved paths, surface worker errors, tick progress
    total_files = 0
    total = len(tasks)
    for completed, (_task, ok, result) in enumerate(run_tasks(
        tasks=tasks, worker=process_single_task, workers=plot_settings.workers,
    ), start=1):
        if not ok:
            on_warning(f"Plot generation failed: {result}")
        elif not result.success:
            on_warning(f"Plot generation failed: {result.error}")
        else:
            for path in (result.png_path, result.svg_path, result.pdf_path, result.html_path):
                if path:
                    saved_files.append(path)
                    total_files += 1
        on_progress(completed, total)
    return total_files


def _run_summary_plots(
    comparison_frame: pd.DataFrame | None,
    root_dir: Path,
    export_formats: list[str],
    saved_files: list[str],
) -> int:
    """Render dual-axis / ridgeline / heatmap summaries from the comparison table.
    Summaries are picture-only — html is filtered out, and an empty format set falls
    back to ["png"]. Appends every written path to saved_files and returns the
    count; returns 0 when there is no comparison frame or it is empty.
    """
    if comparison_frame is None or comparison_frame.empty:
        return 0
    summary_formats = [f for f in export_formats if f in _PICTURE_FORMATS] or ["png"]
    results = generate_all_summary_plots(
        source=comparison_frame,
        output_dir=str(root_dir),
        formats=summary_formats,
    )
    paths = [p for fmt_paths in results.values() for p in fmt_paths.values()]
    saved_files.extend(paths)
    return len(paths)


def run_export(
    *,
    export_directory: str,
    project_name: str,
    run_settings: dict,
    frames: dict[str, Any],
    annotated_queries: AnnotatedQueries,
    reference_store: dict,
    variant_groups: dict,
    comparison_frame: pd.DataFrame | None = None,
    selection: ExportSelection | None = None,
    plot_settings: PlotSettings | None = None,
    callbacks: ExportCallbacks | None = None,
) -> ExportResult:
    """Run a full export end-to-end: versioned dir, parquets, CSV/XLSX, plots, summaries.
    Keyword-only entry point shared by both frontends. Allocates (or reuses) a versioned
    folder for project_name under export_directory, writes run_settings.json,
    then runs each artefact stage gated by selection. run_settings.json is always
    counted as saved even when every stage is disabled.
    Returns an ExportResult whose saved_files lists every written path and whose
    mirror_plot_count / summary_plot_count report per-group totals. Stage failures
    that are recoverable (a bad plot worker, a missing XLSX engine) are surfaced through
    callbacks.on_warning rather than raised.
    """
    # Normalise inputs: defaults for missing dataclasses, no-op for missing callbacks
    selection = selection or ExportSelection()
    plot_settings = plot_settings or PlotSettings()
    callbacks = callbacks or ExportCallbacks()
    on_step = callbacks.on_step or (lambda *_: None)
    on_progress = callbacks.on_plot_progress or (lambda *_: None)
    on_warning = callbacks.on_warning or (lambda *_: None)

    # Stage 1 — Plan & commit folder: reuse on fingerprint match, else allocate next NNN
    root_dir = find_or_create_export_directory(export_directory, project_name, run_settings)
    commit_export_directory(root_dir, run_settings)

    saved_files: list[str] = [str(root_dir / "run_settings.json")]

    # Stage 2 — Parquet frames: write each of the four LazyFrames as {name}.parquet
    if selection.parquet:
        on_step("Saving Parquet DataFrames...")
        for name in _PARQUET_FRAMES:
            frame = frames.get(name)
            if frame is None:
                continue
            out = root_dir / f"{name}.parquet"
            frame.collect().write_parquet(out)
            saved_files.append(str(out))

    # Stage 3 — Comparison table: CSV unconditional; XLSX best-effort (needs openpyxl)
    if selection.csv and comparison_frame is not None and not comparison_frame.empty:
        on_step("Saving comparison summary...")
        csv_path = root_dir / "comparison_summary.csv"
        comparison_frame.to_csv(csv_path, index=False)
        saved_files.append(str(csv_path))
        xlsx_path = root_dir / "comparison_summary.xlsx"
        try:
            comparison_frame.to_excel(xlsx_path, index=False, engine="openpyxl")
            saved_files.append(str(xlsx_path))
        except ImportError:
            on_warning("openpyxl not installed -- skipping XLSX export")
        except Exception as error:
            on_warning(f"XLSX export failed: {error}")

    # Stage 4 — Mirror plots: one per PSM, rendered across the worker pool
    mirror_plot_count = 0
    if selection.mirror_plots or selection.html:
        on_step("Generating mirror plots...")
        # Normalise annotated_queries: dict {psm_id: (query_psm, annotated_peaks_frame)}
        # or list [(psm_id, query_psm, annotated_peaks_frame)]
        if isinstance(annotated_queries, dict):
            queries_list = [
                (psm_id, query_psm, annotated_peaks_frame)
                for psm_id, (query_psm, annotated_peaks_frame) in annotated_queries.items()
            ]
        else:
            queries_list = list(annotated_queries)
        mirror_plot_count = _run_plot_tasks(
            queries_list, reference_store, variant_groups, root_dir,
            plot_settings, selection.html, saved_files, on_progress, on_warning,
        )

    # Stage 5 — Summary plots: dual-axis / ridgeline / heatmap from the comparison frame
    summary_plot_count = 0
    if selection.summary_plots:
        on_step("Generating summary plots...")
        summary_plot_count = _run_summary_plots(
            comparison_frame, root_dir, plot_settings.export_formats, saved_files,
        )

    return ExportResult(
        export_root=root_dir,
        saved_files=saved_files,
        mirror_plot_count=mirror_plot_count,
        summary_plot_count=summary_plot_count,
    )
