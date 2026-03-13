"""
Peptide search functionality for idXML files.

Searches for target peptides and extracts all metadata needed for mzML resolution.
Uses PyOpenMS directly for accurate id_merge_index and spectra_data extraction.
"""
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import numpy as np
import pyopenms as oms

logger = logging.getLogger(__name__)

# Compact abbreviations for common post-translational modifications
MOD_ABBREVIATIONS = {
    "Oxidation": "Ox",
    "Phospho": "Ph",
    "Carbamidomethyl": "CAM",
    "Acetyl": "Ac",
    "Deamidated": "Deam",
    "Methylation": "Me",
    "Dimethyl": "DiMe",
    "Ubiquitin": "GG",
    "Pyro-carbamidomethyl": "Pyro-CAM",
    "Pyro-glu": "Pyro",
}


def abbreviate_mod(name: str) -> str:
    """Return short abbreviation for a modification, or the full name if unknown."""
    return MOD_ABBREVIATIONS.get(name, name)


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
    import importlib
    mod = importlib.import_module("src.peptide_search")
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


def decode_metadata_value(value) -> str:
    """Decode metadata value that might be bytes, list of bytes, or string."""
    if isinstance(value, list) and len(value) > 0:
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value) if value else ""


def extract_spectra_data_map(protein_ids: List) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Extract spectra_data mapping from ProteinIdentification objects.

    Returns:
        Tuple of:
        - run_map: identifier -> list of mzML bare filenames
        - run_map_full: identifier -> list of full metadata paths (as stored in idXML)
    """
    run_map = {}
    run_map_full = {}

    for protein_id in protein_ids:
        identifier = protein_id.getIdentifier()
        mzml_files = []
        mzml_full_paths = []

        # Try spectra_data first (standard for merged idXML)
        if protein_id.metaValueExists("spectra_data"):
            raw_val = protein_id.getMetaValue("spectra_data")
            if isinstance(raw_val, list):
                for v in raw_val:
                    full = decode_metadata_value(v)
                    mzml_full_paths.append(full)
                    mzml_files.append(Path(full).name)
            else:
                full = decode_metadata_value(raw_val)
                mzml_full_paths.append(full)
                mzml_files.append(Path(full).name)

        # Fallback to other metadata keys
        if not mzml_files:
            for key in ["CometAdapter:1:in", "mzML_file", "spectrum_reference"]:
                if protein_id.metaValueExists(key):
                    val = decode_metadata_value(protein_id.getMetaValue(key))
                    if val:
                        mzml_full_paths = [val]
                        mzml_files = [Path(val).name]
                        break

        if mzml_files:
            run_map[identifier] = mzml_files
            run_map_full[identifier] = mzml_full_paths

    return run_map, run_map_full


def get_modification_string(peptide_hit) -> str:
    """
    Extract compact modification string from a PeptideHit.

    Returns format like "M5[Ox]" or "M5[Ox]|S3[Ph]" for multiple mods.
    Terminal mods: "NTerm[Ac]", "CTerm[...]".
    """
    try:
        aa_seq = peptide_hit.getSequence()
        if not aa_seq.isModified():
            return ""

        mods = []
        # Check N-term modification
        if aa_seq.hasNTerminalModification():
            mod = aa_seq.getNTerminalModification()
            if mod:
                mod_name = mod.getFullId().split(" ")[0]  # Strip trailing details
                mods.append(f"NTerm[{abbreviate_mod(mod_name)}]")

        # Check each residue
        for i in range(aa_seq.size()):
            residue = aa_seq.getResidue(i)
            if residue.isModified():
                mod_name = residue.getModificationName()
                abbrev = abbreviate_mod(mod_name)
                mods.append(f"{residue.getOneLetterCode()}{i+1}[{abbrev}]")

        # Check C-term modification
        if aa_seq.hasCTerminalModification():
            mod = aa_seq.getCTerminalModification()
            if mod:
                mod_name = mod.getFullId().split(" ")[0]
                mods.append(f"CTerm[{abbreviate_mod(mod_name)}]")

        return "|".join(mods) if mods else ""
    except Exception as e:
        logger.debug(f"Could not extract modifications: {e}")
        return ""


def _search_single_idxml(
    idxml_path: str,
    idxml_name: str,
    target_set: Set[str],
    mzml_path_strs: List[str],
) -> Tuple[str, List[Tuple[str, PeptideMatch]], List[str]]:
    """
    Worker: search one idXML file for target peptides.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_name, list of (variant_key, PeptideMatch), list of errors)
    """
    from src.path_utils import resolve_mzml_path

    matches: List[Tuple[str, PeptideMatch]] = []
    errors: List[str] = []

    try:
        try:
            from psm_utils.io.idxml import IdXMLReader
            _reader = IdXMLReader(str(idxml_path))
            protein_ids = _reader.protein_ids
            peptide_ids = _reader.peptide_ids
        except Exception as _load_err:
            errors.append(f"No peptide identifications in {idxml_name}: {_load_err}")
            return idxml_name, matches, errors

        if not peptide_ids:
            errors.append(f"No peptide identifications in {idxml_name}")
            return idxml_name, matches, errors

        run_map, run_map_full = extract_spectra_data_map(protein_ids)

        for pep_id in peptide_ids:
            spectrum_id = ""
            if pep_id.metaValueExists("spectrum_reference"):
                spectrum_id = decode_metadata_value(pep_id.getMetaValue("spectrum_reference"))

            id_merge_index = -1
            if pep_id.metaValueExists("id_merge_index"):
                try:
                    id_merge_index = int(pep_id.getMetaValue("id_merge_index"))
                except (ValueError, TypeError):
                    pass

            pep_identifier = pep_id.getIdentifier()
            matched_run_id = None

            if pep_identifier in run_map:
                matched_run_id = pep_identifier
            else:
                for run_id in run_map.keys():
                    if pep_identifier.startswith(run_id) or run_id.startswith(pep_identifier):
                        matched_run_id = run_id
                        break

            mzml_filename = ""
            mzml_metadata_path = ""
            if matched_run_id and matched_run_id in run_map:
                mzml_list = run_map[matched_run_id]
                full_list = run_map_full.get(matched_run_id, [])
                if id_merge_index >= 0 and id_merge_index < len(mzml_list):
                    mzml_filename = mzml_list[id_merge_index]
                    if id_merge_index < len(full_list):
                        mzml_metadata_path = full_list[id_merge_index]
                elif len(mzml_list) == 1:
                    mzml_filename = mzml_list[0]
                    if full_list:
                        mzml_metadata_path = full_list[0]

            hits = pep_id.getHits()
            for hit in hits:
                try:
                    aa_seq = hit.getSequence()
                    base_seq = aa_seq.toUnmodifiedString().upper()

                    if base_seq not in target_set:
                        continue

                    peptidoform_str = aa_seq.toString()
                    charge = hit.getCharge() if hit.getCharge() > 0 else 0
                    score = hit.getScore() if hit.getScore() else 0.0

                    hit_merge_index = id_merge_index
                    if hit.metaValueExists("id_merge_index"):
                        try:
                            hit_merge_index = int(hit.getMetaValue("id_merge_index"))
                        except (ValueError, TypeError):
                            pass

                    hit_mzml_filename = mzml_filename
                    hit_metadata_path = mzml_metadata_path
                    if hit_merge_index != id_merge_index and matched_run_id and matched_run_id in run_map:
                        mzml_list = run_map[matched_run_id]
                        full_list = run_map_full.get(matched_run_id, [])
                        if hit_merge_index >= 0 and hit_merge_index < len(mzml_list):
                            hit_mzml_filename = mzml_list[hit_merge_index]
                            if hit_merge_index < len(full_list):
                                hit_metadata_path = full_list[hit_merge_index]

                    resolved_path = ""
                    if mzml_path_strs and hit_mzml_filename:
                        resolved = resolve_mzml_path(
                            hit_mzml_filename, hit_metadata_path,
                            mzml_path_strs, idxml_path
                        )
                        if resolved:
                            resolved_path = resolved

                    match = PeptideMatch(
                        peptidoform=peptidoform_str,
                        base_sequence=base_seq,
                        spectrum_id=spectrum_id,
                        charge=charge,
                        score=score,
                        idxml_file=idxml_name,
                        mzml_filename=hit_mzml_filename,
                        id_merge_index=hit_merge_index,
                        mzml_resolved_path=resolved_path,
                        idxml_resolved_path=idxml_path
                    )

                    mod_key = get_modification_string(hit)
                    variant_key = f"{base_seq}|{mod_key}" if mod_key else base_seq
                    matches.append((variant_key, match))

                except Exception as e:
                    logger.debug(f"Error processing peptide hit: {e}")
                    continue

    except Exception as e:
        errors.append(f"Error loading {Path(idxml_path).name}: {e}")

    return idxml_name, matches, errors


def _scan_single_idxml_peptides(idxml_path: str) -> Tuple[str, Set[str]]:
    """
    Worker: extract all unique base sequences from one idXML file.

    Must be a module-level function for ProcessPoolExecutor pickling.

    Returns:
        Tuple of (idxml_path, set of base sequences)
    """
    peptides: Set[str] = set()
    try:
        from psm_utils.io.idxml import IdXMLReader
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
    from src.path_utils import disambiguate_paths

    # Normalize target peptides to uppercase, strip whitespace
    target_set = {p.strip().upper() for p in target_peptides if p.strip()}

    if not target_set:
        if error_callback:
            error_callback("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    # Compute disambiguated names for idXML and mzML paths
    idxml_disambig = disambiguate_paths(idxml_paths)
    mzml_path_strs = [str(p) for p in mzml_paths] if mzml_paths else []

    total_files = len(idxml_paths)
    variant_groups: Dict[str, VariantGroup] = {}
    total_psms_found = 0

    def _merge_worker_results(worker_matches, worker_errors):
        """Merge worker results into variant_groups (main process)."""
        nonlocal total_psms_found
        for variant_key, match in worker_matches:
            mod_key = ""
            base_seq = match.base_sequence
            if "|" in variant_key:
                mod_key = variant_key.split("|", 1)[1]

            if variant_key not in variant_groups:
                display_label = format_inline_sequence(base_seq, mod_key)
                variant_groups[variant_key] = VariantGroup(
                    base_sequence=base_seq,
                    modification_key=mod_key,
                    display_label=display_label
                )
            variant_groups[variant_key].matches.append(match)
            total_psms_found += 1

        for err in worker_errors:
            if error_callback:
                error_callback(err)
            logger.error(err)

    workers = max(1, min(workers, total_files))

    if workers <= 1:
        # Sequential path — process files one at a time
        for idx, idxml_path in enumerate(idxml_paths):
            if progress_callback:
                progress_callback(idx, total_files, f"Searching {Path(idxml_path).name}")

            idxml_name = idxml_disambig.get(idxml_path, Path(idxml_path).name)
            _, worker_matches, worker_errors = _search_single_idxml(
                idxml_path, idxml_name, target_set, mzml_path_strs
            )
            _merge_worker_results(worker_matches, worker_errors)

    else:
        # Parallel path — ProcessPoolExecutor with as_completed
        from concurrent.futures import ProcessPoolExecutor, as_completed

        logger.info(f"Parallel idXML search: {total_files} files, {workers} workers")

        with ProcessPoolExecutor(max_workers=workers) as pool:
            future_map = {}
            for idxml_path in idxml_paths:
                idxml_name = idxml_disambig.get(idxml_path, Path(idxml_path).name)
                future = pool.submit(
                    _search_single_idxml,
                    idxml_path, idxml_name, target_set, mzml_path_strs
                )
                future_map[future] = idxml_path

            completed = 0
            for future in as_completed(future_map):
                completed += 1
                idxml_path = future_map[future]
                try:
                    _, worker_matches, worker_errors = future.result()
                    _merge_worker_results(worker_matches, worker_errors)
                except Exception as e:
                    err = f"Worker crashed processing {Path(idxml_path).name}: {e}"
                    if error_callback:
                        error_callback(err)
                    logger.error(err)

                if progress_callback:
                    progress_callback(completed, total_files,
                                      f"Searched {Path(idxml_path).name}")

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
    workers = max(1, min(workers, total_files))

    if workers <= 1:
        # Sequential path
        for idx, idxml_path in enumerate(idxml_paths):
            if progress_callback:
                progress_callback(idx, total_files, f"Scanning {Path(idxml_path).name}")

            if len(unique_peptides) >= limit:
                break

            _, file_peptides = _scan_single_idxml_peptides(idxml_path)
            unique_peptides.update(file_peptides)

    else:
        # Parallel path
        from concurrent.futures import ProcessPoolExecutor, as_completed

        logger.info(f"Parallel peptide scan: {total_files} files, {workers} workers")

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_scan_single_idxml_peptides, path): path
                for path in idxml_paths
            }

            completed = 0
            for future in as_completed(futures):
                completed += 1
                try:
                    _, file_peptides = future.result()
                    unique_peptides.update(file_peptides)
                except Exception as e:
                    logger.warning(f"Worker error scanning peptides: {e}")

                if progress_callback:
                    progress_callback(completed, total_files,
                                      f"Scanned {Path(futures[future]).name}")

                # Cancel remaining futures if limit reached
                if len(unique_peptides) >= limit:
                    for f in futures:
                        f.cancel()
                    break

    if progress_callback:
        progress_callback(total_files, total_files, "Scan complete")

    return sorted(list(unique_peptides))[:limit]


def summarize_search_results(variant_groups: Dict[str, VariantGroup]) -> pd.DataFrame:
    """
    Create a summary DataFrame of search results.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup

    Returns:
        DataFrame with columns: base_sequence, modification, display_label, psm_count, files, mzml_resolved
    """
    rows = []
    for variant_key, vg in variant_groups.items():
        # Get unique files
        idxml_files = set(m.idxml_file for m in vg.matches)
        mzml_files = set(m.mzml_filename for m in vg.matches if m.mzml_filename)
        unresolved = sum(1 for m in vg.matches if not m.mzml_filename)

        rows.append({
            'base_sequence': vg.base_sequence,
            'modification': vg.modification_key or "(none)",
            'display_label': vg.display_label,
            'psm_count': len(vg.matches),
            'idxml_files': len(idxml_files),
            'mzml_resolved': len(vg.matches) - unresolved,
            'unresolved': unresolved
        })

    return pd.DataFrame(rows)
