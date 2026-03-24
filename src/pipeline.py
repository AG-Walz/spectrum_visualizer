"""
Shared processing pipeline helpers for MSpecViz.

Pure functions that are identical between CLI (cli.py) and Streamlit
(content/peptide_mirror.py).  No Streamlit imports — callable from any context.

Extracted functions:
    annotate_spectrum()          – validate + annotate a LoadedSpectrum
    ensure_references_annotated() – annotate a variant's refs once, cache them
    build_comparison_tasks()     – build ComparisonTask list for one query
    write_comparison_summary()   – write CSV / XLSX / summary plots
    mark_skipped_tasks()         – set skip_plot=True for already-exported tasks
    detect_parquet_cache()       – find complete/partial Parquet cache for inputs
    save_parquet_cache()         – save spectra to Parquet with logging
    stream_spectra()             – dispatch loading: Parquet / partial-resume / mzML
    accumulate_plot_result()     – update results_summary + summary_rows from one result
"""
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Spectrum annotation
# ---------------------------------------------------------------------------

def annotate_spectrum(
    spectrum,
    ion_types: List[str],
    mz_tolerance: float,
    add_losses: bool,
    add_first_prefix_ion: bool,
    remove_precursor_peaks: bool,
    normalize: bool,
    errors: List[str],
    log_fn: Optional[Callable[[str], None]] = None,
):
    """
    Validate and annotate a LoadedSpectrum for plotting.

    Returns an annotated DataFrame ready for a ComparisonTask, or None if
    validation fails.  Error messages are appended to the *errors* list and
    (optionally) sent to *log_fn*.

    Args:
        spectrum:              LoadedSpectrum instance
        ion_types:             Ion types to annotate (e.g. ["b", "y"])
        mz_tolerance:          Matching tolerance in Da
        add_losses:            Include neutral loss peaks
        add_first_prefix_ion:  Include first prefix ion
        remove_precursor_peaks: Strip precursor peaks from output
        normalize:             Normalize intensities to [0, 1]
        errors:                Mutable list to which error strings are appended
        log_fn:                Optional callable(str) for additional logging
    """
    from src.peak_annotation import prepare_spectrum_for_plot

    spec_type = "reference" if getattr(spectrum, "is_reference", False) else "query"
    peptidoform = getattr(spectrum, "peptidoform", "<unknown>")

    def _fail(msg: str):
        errors.append(msg)
        if log_fn:
            log_fn(msg)

    if spectrum.mz is None or len(spectrum.mz) == 0:
        _fail(f"Annotation error ({spec_type} {peptidoform}): Empty m/z array")
        return None
    if spectrum.intensity is None or len(spectrum.intensity) == 0:
        _fail(f"Annotation error ({spec_type} {peptidoform}): Empty intensity array")
        return None
    if not peptidoform:
        _fail(f"Annotation error ({spec_type}): Missing peptidoform")
        return None
    if spectrum.precursor_charge <= 0:
        _fail(
            f"Annotation error ({spec_type} {peptidoform}): "
            f"Invalid charge {spectrum.precursor_charge}"
        )
        return None

    try:
        return prepare_spectrum_for_plot(
            mz=spectrum.mz,
            intensity=spectrum.intensity,
            peptidoform=peptidoform,
            charge=spectrum.precursor_charge,
            ion_types=ion_types,
            mz_tolerance=mz_tolerance,
            add_losses=add_losses,
            add_first_prefix_ion=add_first_prefix_ion,
            remove_precursor_peaks=remove_precursor_peaks,
            normalize=normalize,
        )
    except Exception as e:
        _fail(
            f"Annotation error ({spec_type} {peptidoform}, "
            f"charge={spectrum.precursor_charge}, peaks={len(spectrum.mz)}): {e}"
        )
        return None


def ensure_references_annotated(
    variant_key: str,
    refs: list,
    annotated_references: Dict[str, list],
    annotation_kwargs: dict,
    errors: List[str],
    log_fn: Optional[Callable[[str], None]] = None,
) -> list:
    """
    Annotate all reference spectra for a variant, caching the result.

    Modifies *annotated_references* in-place on first call for a variant_key.
    Returns the cached list of (ref, ref_df) pairs.

    Args:
        variant_key:          Peptide variant identifier
        refs:                 List of LoadedSpectrum reference objects
        annotated_references: Shared cache dict (modified in-place)
        annotation_kwargs:    Keyword args forwarded to annotate_spectrum()
        errors:               Mutable list for error accumulation
        log_fn:               Optional logging callable
    """
    if variant_key not in annotated_references:
        annotated_references[variant_key] = []
        for ref in refs:
            ref_df = annotate_spectrum(ref, errors=errors, log_fn=log_fn, **annotation_kwargs)
            if ref_df is not None:
                annotated_references[variant_key].append((ref, ref_df))
    return annotated_references[variant_key]


# ---------------------------------------------------------------------------
# Koina fallback reference
# ---------------------------------------------------------------------------

def _predict_koina_reference(
    query,
    base_sequence: str,
    modification_key: str,
    koina_settings: dict,
) -> Optional["pd.DataFrame"]:
    """
    Generate a Koina-predicted spectrum to use as a mirror-plot reference.

    Returns an annotated DataFrame [mz, intensity, annotation] or None on failure.
    Called from build_comparison_tasks() when no experimental reference exists.
    """
    from src.koina_prediction import predict_koina_spectrum, KOINA_AVAILABLE

    if not KOINA_AVAILABLE:
        return None

    # Reconstruct full peptidoform string for Koina
    peptidoform = getattr(query, "peptidoform", base_sequence)
    charge = query.precursor_charge

    return predict_koina_spectrum(
        peptidoform=str(peptidoform),
        charge=charge,
        model=koina_settings.get("model", "ms2pip_timsTOF2024"),
        collision_energy=koina_settings.get("collision_energy", 35),
        instrument_type=koina_settings.get("instrument_type", "TIMSTOF"),
    )


# ---------------------------------------------------------------------------
# Task building
# ---------------------------------------------------------------------------

def build_comparison_tasks(
    query,
    ref_data: List[Tuple],
    query_df,
    base_sequence: str,
    modification_key: str,
    png_dir: str,
    svg_dir: str,
    pdf_dir: str,
    export_formats: List[str],
    figsize: Tuple,
    dpi: int,
    normalize: bool,
    plot_unpaired: bool,
    task_counter: List[int],
    html_dir: str = "",
    koina_settings: Optional[dict] = None,
) -> list:
    """
    Build ComparisonTask objects for one query spectrum.

    Returns a list of ComparisonTask instances:
      - Empty list if the query has no references, Koina is disabled, and
        plot_unpaired is False
      - One "mirror" task with Koina-predicted reference if no experimental
        reference exists and Koina is enabled
      - One "single" task if no references exist, Koina is disabled, and
        plot_unpaired is True
      - One "mirror" task per reference spectrum otherwise

    Args:
        query:           LoadedSpectrum for the query
        ref_data:        List of (ref_LoadedSpectrum, ref_DataFrame) pairs
        query_df:        Annotated DataFrame for the query
        base_sequence:   Base peptide sequence
        modification_key: Modification description (empty string if none)
        png_dir:         Absolute path for PNG outputs
        svg_dir:         Absolute path for SVG outputs
        pdf_dir:         Absolute path for PDF outputs
        export_formats:  List of format strings, e.g. ["png", "svg"]
        figsize:         (width, height) tuple for matplotlib
        dpi:             Figure DPI
        normalize:       Normalize intensities
        plot_unpaired:   Generate single-spectrum plots when no reference
        task_counter:    Single-element list [int] used as a mutable counter
        html_dir:        Directory for interactive HTML exports
        koina_settings:  Optional dict with keys: enabled, model, collision_energy,
                         instrument_type. When enabled and no experimental reference
                         exists, a Koina-predicted spectrum is used as reference.
    """
    from src.parallel_processor import ComparisonTask
    from src.export import sanitize_filename, shorten_stem
    from src.peptide_search import format_inline_sequence

    tasks: list = []

    if modification_key:
        folder_name = f"{base_sequence}_{sanitize_filename(modification_key)}"
        plot_title = format_inline_sequence(base_sequence, modification_key)
    else:
        folder_name = base_sequence
        plot_title = base_sequence

    if not ref_data:
        # --- Koina fallback: predict a reference spectrum ---
        koina_ref_df = None
        if koina_settings and koina_settings.get("enabled"):
            koina_ref_df = _predict_koina_reference(
                query, base_sequence, modification_key, koina_settings,
            )

        if koina_ref_df is not None:
            # Mirror task with Koina-predicted reference
            query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
            query_label = (
                f"Query: {Path(query.source_file).stem}"
                f"  (RT: {query_rt}, ID: {query.spectrum_id})"
            )
            model_name = koina_settings.get("model", "Koina")
            ref_label = (
                f"Koina: {model_name}"
                f"  (CE={koina_settings.get('collision_energy', 35)}, "
                f"{koina_settings.get('instrument_type', 'TIMSTOF')})"
            )
            query_stem = shorten_stem(Path(query.source_file).stem)
            safe_filename = sanitize_filename(
                f"{base_sequence}_{query_stem}_vs_koina"
            )

            html_plot_path = ""
            if html_dir and "html" in export_formats:
                html_plot_path = str(
                    Path(html_dir) / folder_name / f"{safe_filename}.html"
                )

            tasks.append(ComparisonTask(
                task_id=task_counter[0],
                task_type="mirror",
                query_df=query_df,
                ref_df=koina_ref_df,
                query_label=query_label,
                ref_label=ref_label,
                plot_title=plot_title,
                folder_name=folder_name,
                safe_filename=safe_filename,
                png_dir=png_dir,
                svg_dir=svg_dir,
                pdf_dir=pdf_dir,
                html_plot_path=html_plot_path,
                export_formats=export_formats,
                figsize=figsize,
                dpi=dpi,
                normalize=normalize,
                query_file=Path(query.source_file).name,
                ref_file=f"koina:{model_name}",
                peptide=base_sequence,
                modification=modification_key or "",
                charge=query.precursor_charge,
            ))
            task_counter[0] += 1
            return tasks

        # No Koina reference available — fall back to single-spectrum plot
        if not plot_unpaired:
            return tasks

        # Single-spectrum task (no reference)
        query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
        single_label = (
            f"{Path(query.source_file).stem}"
            f"  (RT: {query_rt}, ID: {query.spectrum_id})"
        )
        query_stem = shorten_stem(Path(query.source_file).stem)
        safe_filename = sanitize_filename(f"{base_sequence}_{query_stem}")

        html_plot_path = ""
        if html_dir and "html" in export_formats:
            html_plot_path = str(
                Path(html_dir) / folder_name / f"{safe_filename}.html"
            )

        tasks.append(ComparisonTask(
            task_id=task_counter[0],
            task_type="single",
            query_df=query_df,
            query_label=single_label,
            plot_title=plot_title,
            folder_name=folder_name,
            safe_filename=safe_filename,
            png_dir=png_dir,
            svg_dir=svg_dir,
            pdf_dir=pdf_dir,
            html_plot_path=html_plot_path,
            export_formats=export_formats,
            figsize=figsize,
            dpi=dpi,
            normalize=normalize,
            query_file=Path(query.source_file).name,
            peptide=base_sequence,
            modification=modification_key or "",
            charge=query.precursor_charge,
        ))
        task_counter[0] += 1
        return tasks

    # Mirror tasks — one per reference
    query_stem = shorten_stem(Path(query.source_file).stem)
    query_rt = f"{query.rt:.1f}" if query.rt else "N/A"
    query_label = (
        f"Query: {Path(query.source_file).stem}"
        f"  (RT: {query_rt}, ID: {query.spectrum_id})"
    )

    for ref, ref_df in ref_data:
        ref_rt = f"{ref.rt:.1f}" if ref.rt else "N/A"
        ref_label = (
            f"Ref: {Path(ref.source_file).stem}"
            f"  (RT: {ref_rt}, ID: {ref.spectrum_id})"
        )
        ref_stem = shorten_stem(Path(ref.source_file).stem)
        safe_filename = sanitize_filename(
            f"{base_sequence}_{query_stem}_vs_{ref_stem}"
        )

        html_plot_path = ""
        if html_dir and "html" in export_formats:
            html_plot_path = str(
                Path(html_dir) / folder_name / f"{safe_filename}.html"
            )

        tasks.append(ComparisonTask(
            task_id=task_counter[0],
            task_type="mirror",
            query_df=query_df,
            ref_df=ref_df,
            query_label=query_label,
            ref_label=ref_label,
            plot_title=plot_title,
            folder_name=folder_name,
            safe_filename=safe_filename,
            png_dir=png_dir,
            svg_dir=svg_dir,
            pdf_dir=pdf_dir,
            html_plot_path=html_plot_path,
            export_formats=export_formats,
            figsize=figsize,
            dpi=dpi,
            normalize=normalize,
            query_file=Path(query.source_file).name,
            ref_file=Path(ref.source_file).name,
            peptide=base_sequence,
            modification=modification_key or "",
            charge=query.precursor_charge,
        ))
        task_counter[0] += 1

    return tasks


# ---------------------------------------------------------------------------
# Plot skip detection
# ---------------------------------------------------------------------------

def mark_skipped_tasks(tasks: list) -> int:
    """
    Mark tasks whose output files already exist on disk as skip_plot=True.

    Must be called before tmpfs path rewriting (iter_chunks) so that paths
    still point to the final disk location.

    Returns the number of tasks skipped.
    """
    skipped = 0
    for task in tasks:
        if task.outputs_exist():
            task.skip_plot = True
            skipped += 1
    return skipped


# ---------------------------------------------------------------------------
# Summary export
# ---------------------------------------------------------------------------

def write_comparison_summary(
    summary_rows: List[dict],
    root_dir: Path,
    export_formats: List[str],
    dpi: int,
    project_name: str,
    log_fn: Optional[Callable[[str], None]] = None,
) -> dict:
    """
    Write comparison summary CSV, XLSX, and summary plots.

    Args:
        summary_rows:    List of per-comparison metric dicts from ComparisonResult
        root_dir:        Export root directory (CSV/XLSX written here)
        export_formats:  Formats for summary plots (e.g. ["png", "svg"])
        dpi:             DPI for summary plots
        project_name:    Used as title prefix for summary plots
        log_fn:          Optional callable(str) for progress messages

    Returns:
        Dict with any of: "csv", "xlsx", "summary_plots"
    """
    import pandas as pd
    from src.spectrum_metrics import shorten_filenames

    if not summary_rows:
        return {}

    result: dict = {}

    try:
        summary_df = pd.DataFrame(summary_rows)

        # Shorten filenames for the display columns
        query_files = summary_df["query_file"].tolist() if "query_file" in summary_df.columns else []
        ref_files   = summary_df["ref_file"].tolist()   if "ref_file"   in summary_df.columns else []
        all_filenames = list(set(query_files + ref_files))
        short_map = dict(zip(all_filenames, shorten_filenames(all_filenames)))

        if "query_file" in summary_df.columns:
            summary_df.insert(0, "query_short", summary_df["query_file"].map(short_map))
        if "ref_file" in summary_df.columns:
            summary_df.insert(1, "ref_short", summary_df["ref_file"].map(short_map))

        # Core columns + any ims_ supplementary stats
        from src.parallel_processor import _IMS_STAT_KEYS
        ims_cols = [f"ims_{k}" for k in _IMS_STAT_KEYS]

        col_order = [
            "query_short", "ref_short", "peptide", "modification", "charge",
            "R2m", "R2b", "Cos", "Ent", "Dice", "Jacc", "ShR",
            "total_q", "total_r", "ann_q", "ann_r", "shared",
            *ims_cols,
            "query_file", "ref_file",
        ]
        summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]

        rename_map = {
            "query_short": "Query",       "ref_short": "Reference",
            "peptide": "Peptide",         "modification": "Modification",
            "charge": "Charge",
            "R2m": "R² match",            "R2b": "R² binned",
            "Cos": "Cosine",              "Ent": "Entropy",
            "Dice": "Dice",               "Jacc": "Jaccard",   "ShR": "SharedRatio",
            "total_q": "Peaks (Q)",       "total_r": "Peaks (R)",
            "ann_q": "Annotated (Q)",     "ann_r": "Annotated (R)",
            "shared": "Shared Ions",
            "query_file": "Query File",   "ref_file": "Ref File",
        }
        summary_df = summary_df.rename(columns=rename_map)

        # CSV (primary)
        csv_path = root_dir / "comparison_summary.csv"
        summary_df.to_csv(csv_path, index=False)
        result["csv"] = str(csv_path)
        if log_fn:
            log_fn(f"Summary CSV: {csv_path}")

        # XLSX (secondary — optional openpyxl dependency)
        try:
            xlsx_path = root_dir / "comparison_summary.xlsx"
            summary_df.to_excel(xlsx_path, index=False, sheet_name="Comparisons")
            result["xlsx"] = str(xlsx_path)
            if log_fn:
                log_fn(f"Summary XLSX: {xlsx_path}")
        except Exception as xlsx_err:
            if log_fn:
                log_fn(f"Warning: Could not write XLSX: {xlsx_err}")

        # Summary plots (optional)
        try:
            from src.summary_plots import generate_all_summary_plots
            summary_plot_paths = generate_all_summary_plots(
                source=summary_df,
                output_dir=str(root_dir),
                metric="R² binned",
                formats=export_formats,
                dpi=dpi,
                title_prefix=project_name,
            )
            result["summary_plots"] = summary_plot_paths
            if log_fn:
                n_sp = sum(len(v) for v in summary_plot_paths.values())
                log_fn(f"Summary plots: {n_sp} files")
        except Exception as sp_err:
            if log_fn:
                log_fn(f"Warning: Could not generate summary plots: {sp_err}")

    except Exception as e:
        if log_fn:
            log_fn(f"Warning: Could not write summary: {e}")
        logger.warning("Summary export failed", exc_info=True)

    return result


# ---------------------------------------------------------------------------
# Parquet cache detection
# ---------------------------------------------------------------------------

@dataclass
class CacheDetectionResult:
    """Result of auto-detecting an existing Parquet cache."""
    cache_fingerprint: Optional[str]   # None if detection was skipped
    cache_dir: Optional[str]           # None if detection was skipped
    parquet_path: Optional[str]        # Set when a complete cache was found
    partial_cache_path: Optional[str]  # Set when only an incomplete cache exists
    loaded_from_cache: bool            # True when parquet_path was auto-detected


def detect_parquet_cache(
    idxml_files: List[str],
    mzml_files: List[str],
    peptides: List[str],
    reference_substring: str,
    export_dir: str,
    cache_subdir: str,
    log_fn: Optional[Callable[[str], None]] = None,
    idxml_filter: str = "",
    mzml_filter: str = "",
) -> CacheDetectionResult:
    """
    Detect an existing Parquet cache for this exact set of inputs.

    Computes a fingerprint from (idxml_files, mzml_files, peptides,
    reference_substring, idxml_filter, mzml_filter) and scans *cache_dir*
    for a matching file.

    Returns a CacheDetectionResult with:
      - ``parquet_path``       set if a complete cache was found → skip mzML loading
      - ``partial_cache_path`` set if an incomplete cache exists → resume from it
      - ``loaded_from_cache``  True when parquet_path was auto-detected

    Args:
        idxml_files:        List of idXML file paths (inputs to fingerprint)
        mzml_files:         List of mzML file paths (inputs to fingerprint)
        peptides:           Target peptide list (inputs to fingerprint)
        reference_substring: Substring identifying reference files
        export_dir:         Root export directory
        cache_subdir:       Sub-directory name for cache files (e.g. ".mspecviz_cache")
        log_fn:             Optional callable(str) for progress/status messages
        idxml_filter:       idXML filter string (included in fingerprint)
        mzml_filter:        mzML filter string (included in fingerprint)
    """
    from src.parquet_cache import (
        compute_cache_fingerprint,
        find_existing_cache,
        get_parquet_row_counts,
    )

    cache_dir = str(Path(export_dir) / cache_subdir)
    cache_fingerprint = compute_cache_fingerprint(
        idxml_files, mzml_files, peptides, reference_substring,
        idxml_filter=idxml_filter, mzml_filter=mzml_filter,
    )

    # Complete cache?
    cached_path = find_existing_cache(cache_dir, cache_fingerprint)
    if cached_path:
        if log_fn:
            log_fn(f"Cache found! Loading from: {Path(cached_path).name}")
        return CacheDetectionResult(
            cache_fingerprint=cache_fingerprint,
            cache_dir=cache_dir,
            parquet_path=cached_path,
            partial_cache_path=None,
            loaded_from_cache=True,
        )

    # Partial (incomplete) cache?
    partial_path = find_existing_cache(cache_dir, cache_fingerprint, allow_partial=True)
    if partial_path:
        p_refs, p_queries = get_parquet_row_counts(partial_path)
        if log_fn:
            log_fn(
                f"Partial cache found ({p_refs} refs, {p_queries} queries). "
                f"Resuming — loading remaining spectra from mzML."
            )
        return CacheDetectionResult(
            cache_fingerprint=cache_fingerprint,
            cache_dir=cache_dir,
            parquet_path=None,
            partial_cache_path=partial_path,
            loaded_from_cache=False,
        )

    return CacheDetectionResult(
        cache_fingerprint=cache_fingerprint,
        cache_dir=cache_dir,
        parquet_path=None,
        partial_cache_path=None,
        loaded_from_cache=False,
    )


# ---------------------------------------------------------------------------
# Parquet cache save
# ---------------------------------------------------------------------------

def save_parquet_cache(
    cached_references: Dict[str, list],
    cached_queries: list,
    cache_path: str,
    cache_metadata: dict,
    cache_complete: List[bool],
    partial: bool = False,
    log_fn: Optional[Callable[[str], None]] = None,
) -> bool:
    """
    Save accumulated spectra to a Parquet cache file.

    Handles the save call, success/error logging, and marks *cache_complete[0]*
    True on a successful complete (non-partial) save.

    Args:
        cached_references: Dict mapping variant_key → list of LoadedSpectrum
        cached_queries:    List of (variant_key, LoadedSpectrum) pairs
        cache_path:        Absolute path for the output .parquet file
        cache_metadata:    Metadata dict written into the Parquet file
        cache_complete:    Single-element list [bool] — set to True on success
        partial:           If True, marks the cache as incomplete
        log_fn:            Optional callable(str) for progress messages

    Returns:
        True if saved successfully, False otherwise.
    """
    from src.parquet_cache import save_spectra_to_parquet

    try:
        saved = save_spectra_to_parquet(
            reference_spectra=cached_references,
            query_spectra=cached_queries,
            output_path=cache_path,
            cache_metadata=cache_metadata,
            partial=partial,
        )
        if saved:
            if not partial:
                cache_complete[0] = True
            n_refs = sum(len(v) for v in cached_references.values())
            tag = "partial" if partial else "complete"
            if log_fn:
                log_fn(
                    f"Cache saved [{tag}]: {Path(saved).name} "
                    f"({len(cached_queries)} queries, {n_refs} refs)"
                )
            return True
    except Exception as cache_err:
        if log_fn:
            log_fn(f"Warning: cache save failed: {cache_err}")
        logger.warning("Cache save failed", exc_info=True)
    return False


# ---------------------------------------------------------------------------
# Spectrum loading dispatch
# ---------------------------------------------------------------------------

def stream_spectra(
    use_parquet: bool,
    parquet_file: Optional[str],
    partial_cache_path: Optional[str],
    variant_groups: Optional[dict],
    mzml_paths: Optional[list],
    reference_substring: str,
    workers: int,
    on_reference_loaded: Callable,
    on_query_ready: Callable,
    on_file_start: Optional[Callable],
    on_file_complete: Optional[Callable],
    on_error: Callable,
    cached_references: Dict[str, list],
    log_fn: Optional[Callable[[str], None]] = None,
) -> None:
    """
    Dispatch spectrum loading to the appropriate backend and populate cached_references.

    Three paths — selected by the combination of use_parquet / partial_cache_path:

    1. **Parquet fast path** (use_parquet=True):
       Load all spectra from *parquet_file* via callbacks.

    2. **Partial-resume path** (partial_cache_path set):
       First replay the partial Parquet through callbacks, then stream the
       remaining spectra from mzML files (skipping already-cached IDs).

    3. **Fresh mzML path** (neither):
       Stream all spectra from mzML files.

    *cached_references* is updated in-place with any reference spectra encountered.

    Args:
        use_parquet:          True when loading from a complete Parquet file
        parquet_file:         Path to the Parquet file (required when use_parquet)
        partial_cache_path:   Path to an incomplete Parquet file (resume path)
        variant_groups:       Variant → VariantGroup mapping (mzML path only)
        mzml_paths:           List of mzML Path objects (mzML path only)
        reference_substring:  Substring identifying reference files
        workers:              Number of parallel workers for mzML streaming
        on_reference_loaded:  Callback(variant_key, count) for loaded references
        on_query_ready:       Callback(StreamingResult) for each ready query
        on_file_start:        Optional callback(filename, is_reference)
        on_file_complete:     Optional callback(filename, loaded, not_found)
        on_error:             Callback(error_msg) for errors
        cached_references:    Mutable dict updated with reference spectra in-place
        log_fn:               Optional callable(str) for status messages
    """
    from src.spectrum_loader import load_spectra_from_parquet, load_spectra_streaming

    if use_parquet:
        if log_fn:
            log_fn(f"Loading from Parquet: {Path(parquet_file).name}")
        refs = load_spectra_from_parquet(
            parquet_path=parquet_file,
            reference_substring=reference_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            error_callback=on_error,
        )
        cached_references.update(refs)

    elif partial_cache_path:
        if log_fn:
            log_fn(f"Resuming from partial cache: {Path(partial_cache_path).name}")

        # Step 1: Replay partial Parquet through callbacks
        refs = load_spectra_from_parquet(
            parquet_path=partial_cache_path,
            reference_substring=reference_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            error_callback=on_error,
        )
        cached_references.update(refs)

        # Step 2: Build skip set from partial cache
        from src.parquet_cache import get_cached_spectrum_ids
        skip_ids = get_cached_spectrum_ids(partial_cache_path)
        if log_fn:
            log_fn(f"Skip set: {len(skip_ids)} spectra already cached")

        # Step 3: Stream remaining spectra, skipping already-cached IDs
        new_refs = load_spectra_streaming(
            variant_groups=variant_groups,
            mzml_paths=mzml_paths,
            reference_substring=reference_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=on_file_complete,
            error_callback=on_error,
            workers=workers,
            skip_ids=skip_ids,
        )
        for vk, specs in new_refs.items():
            if specs:
                if vk not in cached_references:
                    cached_references[vk] = []
                cached_references[vk].extend(specs)

    else:
        load_spectra_streaming(
            variant_groups=variant_groups,
            mzml_paths=mzml_paths,
            reference_substring=reference_substring,
            on_reference_loaded=on_reference_loaded,
            on_query_ready=on_query_ready,
            on_file_start=on_file_start,
            on_file_complete=on_file_complete,
            error_callback=on_error,
            workers=workers,
        )


# ---------------------------------------------------------------------------
# Plot result accumulation
# ---------------------------------------------------------------------------

def accumulate_plot_result(
    plot_result,
    results_summary: dict,
    summary_rows: list,
) -> None:
    """
    Update results_summary and summary_rows from a single ComparisonResult.

    Intended to be called from an on_plot_complete callback in both CLI and
    Streamlit, keeping the data-accumulation logic in one place while each
    caller handles its own progress-bar UI.

    Args:
        plot_result:    ComparisonResult from run_parallel_comparisons
        results_summary: Mutable dict with "plots_generated", "files_exported", "errors"
        summary_rows:   Mutable list of metric dicts for summary CSV/XLSX
    """
    if plot_result.success:
        results_summary["plots_generated"] += 1
        if plot_result.png_path:
            results_summary["files_exported"]["png"].append(plot_result.png_path)
        if plot_result.svg_path:
            results_summary["files_exported"]["svg"].append(plot_result.svg_path)
        if plot_result.summary_row:
            summary_rows.append(plot_result.summary_row)
    else:
        results_summary["errors"].append(f"Plot error: {plot_result.error}")
