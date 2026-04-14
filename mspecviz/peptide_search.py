"""
Peptide search functionality for identification files.

Searches for target peptides using psm_utils PSM iteration; psm_utils resolves
id_merge_index and spectra_data internally so the per-PSM loop is format-agnostic.
Parallel/sequential dispatch is delegated to ``mspecviz.common.parallel.run_tasks``;
pyopenms/psm_utils metadata helpers live in ``mspecviz.idxml_parsing``.
"""
import importlib
import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import polars as pl
from psm_utils.io.idxml import IdXMLReader

from mspecviz.common.parallel import run_tasks
from mspecviz.idxml_parsing import (
    MOD_ABBREVIATIONS,
    abbreviate_mod,
    extract_spectra_data_map,
    mod_key_from_peptidoform,
)
from mspecviz.path_utils import disambiguate_paths, resolve_mzml_path

logger = logging.getLogger(__name__)

# Re-exports for backwards-compat with callers that import from peptide_search
__all__ = [
    "MOD_ABBREVIATIONS", "abbreviate_mod", "format_inline_sequence",
    "PeptideMatch", "VariantGroup",
    "search_peptides_in_idxmls", "get_unique_peptides_from_idxmls",
    "count_psms_per_peptide", "summarize_search_results",
]


def format_inline_sequence(base_sequence: str, modification_key: str) -> str:
    """
    Build an inline-annotated sequence from base_sequence and compact modification_key.

    Args:
        base_sequence: Unmodified sequence, e.g. "PEPTMDEK"
        modification_key: Compact mod key, e.g. "M5[Ox]" or "M5[Ox]|S3[Ph]"
                          Terminal mods: "NTerm[Ac]", "CTerm[...]"

    Returns:
        Inline string, e.g. "PEPTM[Ox]DEK" or "PEPTMDEK" if no mods
    """
    if not modification_key:
        return base_sequence

    # Parse each mod entry: "M5[Ox]" -> (position=5, abbrev="Ox")
    # Terminal mods: "NTerm[Ac]" -> (position=0, abbrev="Ac", is_nterm=True)
    insertions = {}  # position (1-based) -> list of abbreviations to insert after that position
    prefix_mods = []  # N-terminal mods shown as prefix
    suffix_mods = []  # C-terminal mods shown as suffix

    for entry in modification_key.split("|"):
        entry = entry.strip()
        if not entry:
            continue

        # Extract abbreviation from brackets
        bracket_start = entry.find("[")
        bracket_end = entry.find("]")
        if bracket_start < 0 or bracket_end < 0:
            continue

        abbrev = entry[bracket_start:bracket_end + 1]  # e.g. "[Ox]"
        prefix = entry[:bracket_start]  # e.g. "M5" or "NTerm" or "CTerm"

        if prefix == "NTerm":
            prefix_mods.append(abbrev)
        elif prefix == "CTerm":
            suffix_mods.append(abbrev)
        else:
            # Parse residue + position: "M5" -> position 5
            # Find where digits start
            pos_str = ""
            for ch in prefix:
                if ch.isdigit():
                    pos_str += ch
            if pos_str:
                pos = int(pos_str)
                if pos not in insertions:
                    insertions[pos] = []
                insertions[pos].append(abbrev)

    # Build the inline sequence
    parts = []
    for part in prefix_mods:
        parts.append(part)
    for i, aa in enumerate(base_sequence):
        parts.append(aa)
        pos = i + 1  # 1-based position
        if pos in insertions:
            for abbrev in insertions[pos]:
                parts.append(abbrev)
    for part in suffix_mods:
        parts.append(part)

    return "".join(parts)


def _rebuild_peptide_match(*args):
    """Reconstruct PeptideMatch from pickled args (avoids Streamlit reload identity issues)."""
    mod = importlib.import_module("mspecviz.peptide_search")
    return mod.PeptideMatch(*args)


@dataclass
class PeptideMatch:
    """A single peptide-spectrum match with all metadata for mzML resolution."""
    peptidoform: str          # Full peptidoform string with modifications
    base_sequence: str        # Unmodified sequence
    spectrum_id: str          # Native ID in mzML (e.g., "index=12345")
    charge: int               # Precursor charge
    score: float              # PSM score
    idxml_file: str           # Source idXML display name (disambiguated)
    # Direct mzML mapping from PyOpenMS metadata
    mzml_filename: str        # Bare mzML filename from id_merge_index
    id_merge_index: int       # Index into spectra_data list (-1 if unavailable)
    # Resolved full paths for disambiguation (set when mzml_paths are provided)
    mzml_resolved_path: str = ""   # Full resolved path to the mzML file
    idxml_resolved_path: str = ""  # Full path to the source idXML file

    def __reduce__(self):
        """Custom pickle: always resolve via canonical module path."""
        return (_rebuild_peptide_match, (
            self.peptidoform, self.base_sequence, self.spectrum_id,
            self.charge, self.score, self.idxml_file,
            self.mzml_filename, self.id_merge_index,
            self.mzml_resolved_path, self.idxml_resolved_path,
        ))


@dataclass
class VariantGroup:
    """Group of spectra for one peptide variant (sequence + modification)."""
    base_sequence: str
    modification_key: str     # e.g., "M5[Ox]" or "" for unmodified
    display_label: str        # Human-readable label
    matches: List[PeptideMatch] = field(default_factory=list)
    reference_spectra: List = field(default_factory=list)
    query_spectra: List = field(default_factory=list)


def _build_run_to_metadata_map(protein_ids) -> Dict[str, str]:
    """Return {run_stem: full_path} for mzML disambiguation."""
    _, run_map_full = extract_spectra_data_map(protein_ids)
    run_to_metadata: Dict[str, str] = {}
    for full_paths in run_map_full.values():
        for full_path in full_paths:
            stem = Path(full_path).stem
            if stem not in run_to_metadata:
                run_to_metadata[stem] = full_path
    return run_to_metadata


def _build_match_from_psm(
    psm,
    idxml_path: str,
    idxml_name: str,
    run_to_metadata: Dict[str, str],
    mzml_path_strs: List[str],
) -> Optional[Tuple[str, PeptideMatch]]:
    """Build a (variant_key, PeptideMatch) from a psm_utils PSM, or None."""
    base_seq = psm.peptidoform.sequence.upper()

    mzml_filename = (psm.run + ".mzML") if psm.run else ""
    metadata_path = run_to_metadata.get(psm.run, "") if psm.run else ""

    resolved_path = ""
    if mzml_path_strs and mzml_filename:
        resolved = resolve_mzml_path(
            mzml_filename, metadata_path, mzml_path_strs, idxml_path,
        )
        if resolved:
            resolved_path = resolved

    mod_key = mod_key_from_peptidoform(psm.peptidoform)
    variant_key = f"{base_seq}|{mod_key}" if mod_key else base_seq

    match = PeptideMatch(
        peptidoform=psm.peptidoform.proforma.split("/")[0],
        base_sequence=base_seq,
        spectrum_id=psm.spectrum_id or "",
        charge=psm.get_precursor_charge() or 0,
        score=psm.score or 0.0,
        idxml_file=idxml_name,
        mzml_filename=mzml_filename,
        id_merge_index=-1,
        mzml_resolved_path=resolved_path,
        idxml_resolved_path=idxml_path,
    )
    return variant_key, match


def _search_single_idxml(
    idxml_path: str,
    idxml_name: str,
    target_set: Set[str],
    mzml_path_strs: List[str],
) -> Tuple[str, List[Tuple[str, PeptideMatch]], List[str]]:
    """
    Worker: search one idXML file for target peptides using psm_utils PSM iteration.

    psm_utils resolves id_merge_index → run name internally, so each psm.run
    already contains the correct mzML basename (without extension).

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_name, list of (variant_key, PeptideMatch), list of errors)
    """
    matches: List[Tuple[str, PeptideMatch]] = []
    errors: List[str] = []

    try:
        try:
            _reader = IdXMLReader(str(idxml_path))
        except Exception as _load_err:
            errors.append(f"No peptide identifications in {idxml_name}: {_load_err}")
            return idxml_name, matches, errors

        run_to_metadata = _build_run_to_metadata_map(_reader.protein_ids)

        # Manual iteration so a parse error inside psm_utils __iter__ (e.g. float(None) on a
        # missing rescoring feature) is reported clearly. psm_utils uses a generator — once
        # it raises, iteration stops. Affected files are typically unfiltered Percolator
        # output (_pout.idXML).
        psm_iter = iter(_reader)
        psms_read = 0
        parse_error: Optional[str] = None
        while True:
            try:
                psm = next(psm_iter)
                psms_read += 1
            except StopIteration:
                break
            except Exception as e:
                parse_error = str(e)
                logger.debug(f"PSM parse error in {idxml_name} after {psms_read} PSMs: {e}")
                break  # generator is dead — no point continuing

            try:
                if psm.peptidoform.sequence.upper() not in target_set:
                    continue
                built = _build_match_from_psm(
                    psm, idxml_path, idxml_name, run_to_metadata, mzml_path_strs,
                )
                if built is not None:
                    matches.append(built)
            except Exception as e:
                logger.debug(f"Error processing PSM: {e}")
                continue

        if parse_error is not None:
            errors.append(
                f"{idxml_name}: stopped after {psms_read} PSMs — a PSM could not be "
                f"parsed (null rescoring feature value). This happens with unfiltered "
                f"Percolator output (_pout.idXML) that contains PSMs missing metadata. "
                f"Switch to _pout_filtered.idXML for complete, reliable results."
            )

    except Exception as e:
        errors.append(f"Error loading {Path(idxml_path).name}: {e}")

    return idxml_name, matches, errors


def count_psms_per_peptide(idxml_path: str) -> Dict[str, int]:
    """Count PSMs per base peptide sequence in one idXML file.

    Module-level so it can be submitted to a ProcessPoolExecutor.
    Returns an empty dict on any parsing error.
    """
    counts: Counter = Counter()
    try:
        for psm in IdXMLReader(str(idxml_path)):
            seq = psm.peptidoform.sequence.upper()
            if seq:
                counts[seq] += 1
    except Exception:
        pass
    return dict(counts)


def _scan_single_idxml_peptides(idxml_path: str) -> Tuple[str, Set[str]]:
    """
    Worker: extract all unique base sequences from one idXML file.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_path, set of base sequences)
    """
    peptides: Set[str] = set()
    try:
        _reader = IdXMLReader(str(idxml_path))
        for _psm in _reader:
            _base_seq = _psm.peptidoform.sequence.upper()
            if _base_seq:
                peptides.add(_base_seq)
    except Exception:
        pass  # If psm_utils fails, return empty (same as current error handling)

    return idxml_path, peptides


def search_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[callable] = None,
    error_callback: Optional[callable] = None,
    workers: int = 1,
) -> Dict[str, VariantGroup]:
    """
    Search for target peptides across multiple idXML files using PyOpenMS.

    Extracts all metadata needed for mzML resolution including id_merge_index.

    Args:
        idxml_paths: List of idXML file paths to search
        target_peptides: List of peptide sequences to find (base sequences)
        mzml_paths: Optional list of selected mzML file paths (for disambiguation)
        progress_callback: Optional callback(current, total, message) for progress
        error_callback: Optional callback(error_message) for immediate error reporting
        workers: Number of parallel workers. 1 = sequential (default).

    Returns:
        Dict mapping variant_key -> VariantGroup
    """
    target_set = {p.strip().upper() for p in target_peptides if p.strip()}

    if not target_set:
        if error_callback:
            error_callback("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    idxml_disambig = disambiguate_paths(idxml_paths)
    mzml_path_strs = [str(p) for p in mzml_paths] if mzml_paths else []

    total_files = len(idxml_paths)
    variant_groups: Dict[str, VariantGroup] = {}
    total_psms_found = 0

    tasks = [
        (p, idxml_disambig.get(p, Path(p).name), target_set, mzml_path_strs)
        for p in idxml_paths
    ]

    if workers > 1:
        logger.info(f"Parallel idXML search: {total_files} files, {workers} workers")

    def _on_result(task, result):
        nonlocal total_psms_found
        _, worker_matches, worker_errors = result
        for variant_key, match in worker_matches:
            mod_key = variant_key.split("|", 1)[1] if "|" in variant_key else ""
            if variant_key not in variant_groups:
                variant_groups[variant_key] = VariantGroup(
                    base_sequence=match.base_sequence,
                    modification_key=mod_key,
                    display_label=format_inline_sequence(match.base_sequence, mod_key),
                )
            variant_groups[variant_key].matches.append(match)
            total_psms_found += 1
        for err in worker_errors:
            if error_callback:
                error_callback(err)
            logger.error(err)

    def _on_error(task, exc):
        err = f"Worker crashed processing {Path(task[0]).name}: {exc}"
        if error_callback:
            error_callback(err)
        logger.error(err)

    def _on_progress(completed, total):
        if progress_callback:
            progress_callback(completed, total, "Searching idXML files")

    run_tasks(
        tasks=tasks,
        worker=_search_single_idxml,
        workers=workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
    )

    if progress_callback:
        progress_callback(total_files, total_files, "Search complete")

    logger.info(f"Found {total_psms_found} PSMs across {len(variant_groups)} variants")

    return variant_groups


def get_unique_peptides_from_idxmls(
    idxml_paths: List[str],
    limit: int = 100,
    progress_callback: Optional[callable] = None,
    workers: int = 1,
) -> List[str]:
    """
    Extract unique peptide sequences from idXML files.

    Args:
        idxml_paths: List of idXML file paths
        limit: Maximum number of unique peptides to return
        progress_callback: Optional callback for progress
        workers: Number of parallel workers. 1 = sequential (default).

    Returns:
        List of unique base sequences (sorted alphabetically)
    """
    unique_peptides: Set[str] = set()
    total_files = len(idxml_paths)

    if workers > 1:
        logger.info(f"Parallel peptide scan: {total_files} files, {workers} workers")

    def _on_result(task, result):
        _, file_peptides = result
        unique_peptides.update(file_peptides)

    def _on_error(task, exc):
        logger.warning(f"Worker error scanning peptides: {exc}")

    def _on_progress(completed, total):
        if progress_callback:
            progress_callback(completed, total, "Scanning idXML files")

    run_tasks(
        tasks=idxml_paths,  # single-arg → auto-wrapped by run_tasks
        worker=_scan_single_idxml_peptides,
        workers=workers,
        on_result=_on_result,
        on_error=_on_error,
        on_progress=_on_progress,
        stop_when=lambda: len(unique_peptides) >= limit,
    )

    if progress_callback:
        progress_callback(total_files, total_files, "Scan complete")

    return sorted(list(unique_peptides))[:limit]


def summarize_search_results(variant_groups: Dict[str, VariantGroup]) -> pl.LazyFrame:
    """
    Create a variant-level summary LazyFrame for the openms_insight Table component.

    One row per unique peptide variant (base_sequence + modification combination).

    Args:
        variant_groups: Dict of variant_key -> VariantGroup

    Returns:
        Polars LazyFrame with columns:
        - id (Int32): Row identifier required by openms_insight Table
        - base_sequence (String): Unmodified amino acid sequence
        - modification (String): Modification key or "(unmodified)"
        - display_label (String): Inline-annotated sequence (e.g. "PEPTM[Ox]DEK")
        - psm_count (Int32): Number of PSMs found for this variant
        - idxml_files (Int32): Number of distinct idXML files with hits
        - mzml_resolved (Int32): PSMs with a resolved mzML path
        - unresolved (Int32): PSMs where mzML could not be resolved
    """
    rows = []
    for row_id, (variant_key, variant_group) in enumerate(variant_groups.items()):
        idxml_files = set(m.idxml_file for m in variant_group.matches)
        unresolved = sum(1 for m in variant_group.matches if not m.mzml_filename)
        rows.append({
            "id": row_id,
            "base_sequence": variant_group.base_sequence,
            "modification": variant_group.modification_key or "(unmodified)",
            "display_label": variant_group.display_label,
            "psm_count": len(variant_group.matches),
            "idxml_files": len(idxml_files),
            "mzml_resolved": len(variant_group.matches) - unresolved,
            "unresolved": unresolved,
        })

    if not rows:
        schema = {
            "id": pl.Int32, "base_sequence": pl.String, "modification": pl.String,
            "display_label": pl.String, "psm_count": pl.Int32,
            "idxml_files": pl.Int32, "mzml_resolved": pl.Int32, "unresolved": pl.Int32,
        }
        return pl.DataFrame(schema=schema).lazy()

    return pl.DataFrame(rows).with_columns([
        pl.col("id").cast(pl.Int32),
        pl.col("psm_count").cast(pl.Int32),
        pl.col("idxml_files").cast(pl.Int32),
        pl.col("mzml_resolved").cast(pl.Int32),
        pl.col("unresolved").cast(pl.Int32),
    ]).lazy()
