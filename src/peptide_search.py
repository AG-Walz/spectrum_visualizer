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


@dataclass
class VariantGroup:
    """Group of spectra for one peptide variant (sequence + modification)."""
    base_sequence: str
    modification_key: str     # e.g., "Oxidation@M5" or "" for unmodified
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
    """Extract modification string from a PeptideHit."""
    try:
        aa_seq = peptide_hit.getSequence()
        if not aa_seq.isModified():
            return ""

        mods = []
        # Check N-term modification
        if aa_seq.hasNTerminalModification():
            mod = aa_seq.getNTerminalModification()
            if mod:
                mods.append(f"NTerm:{mod.getFullId()}")

        # Check each residue
        for i in range(aa_seq.size()):
            residue = aa_seq.getResidue(i)
            if residue.isModified():
                mod_name = residue.getModificationName()
                mods.append(f"{mod_name}@{residue.getOneLetterCode()}{i+1}")

        # Check C-term modification
        if aa_seq.hasCTerminalModification():
            mod = aa_seq.getCTerminalModification()
            if mod:
                mods.append(f"CTerm:{mod.getFullId()}")

        return "|".join(mods) if mods else ""
    except Exception as e:
        logger.debug(f"Could not extract modifications: {e}")
        return ""


def search_peptides_in_idxmls(
    idxml_paths: List[str],
    target_peptides: List[str],
    mzml_paths: Optional[List[str]] = None,
    progress_callback: Optional[callable] = None,
    error_callback: Optional[callable] = None
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

    Returns:
        Dict mapping variant_key -> VariantGroup
    """
    from src.path_utils import disambiguate_paths, resolve_mzml_path

    # Normalize target peptides to uppercase, strip whitespace
    target_set = {p.strip().upper() for p in target_peptides if p.strip()}

    if not target_set:
        if error_callback:
            error_callback("No valid target peptides provided")
        logger.warning("No valid target peptides provided")
        return {}

    # Compute disambiguated names for idXML and mzML paths
    idxml_disambig = disambiguate_paths(idxml_paths)
    mzml_disambig = disambiguate_paths(mzml_paths) if mzml_paths else {}
    mzml_path_strs = [str(p) for p in mzml_paths] if mzml_paths else []

    variant_groups: Dict[str, VariantGroup] = {}
    total_files = len(idxml_paths)
    total_psms_found = 0

    for idx, idxml_path in enumerate(idxml_paths):
        if progress_callback:
            progress_callback(idx, total_files, f"Searching {Path(idxml_path).name}")

        try:
            # Load idXML with PyOpenMS
            protein_ids = []
            peptide_ids = []
            oms.IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)

            idxml_name = idxml_disambig.get(idxml_path, Path(idxml_path).name)

            if not peptide_ids:
                if error_callback:
                    error_callback(f"No peptide identifications in {idxml_name}")
                continue

            # Extract spectra_data mapping for this file
            run_map, run_map_full = extract_spectra_data_map(protein_ids)

            # Process each PeptideIdentification
            for pep_id in peptide_ids:
                # Get spectrum reference (native ID)
                spectrum_id = ""
                if pep_id.metaValueExists("spectrum_reference"):
                    spectrum_id = decode_metadata_value(pep_id.getMetaValue("spectrum_reference"))

                # Get id_merge_index if available
                id_merge_index = -1
                if pep_id.metaValueExists("id_merge_index"):
                    try:
                        id_merge_index = int(pep_id.getMetaValue("id_merge_index"))
                    except (ValueError, TypeError):
                        pass

                # Find matching run identifier
                pep_identifier = pep_id.getIdentifier()
                matched_run_id = None

                if pep_identifier in run_map:
                    matched_run_id = pep_identifier
                else:
                    # Try prefix matching
                    for run_id in run_map.keys():
                        if pep_identifier.startswith(run_id) or run_id.startswith(pep_identifier):
                            matched_run_id = run_id
                            break

                # Resolve mzML filename
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

                # Process each PeptideHit
                hits = pep_id.getHits()
                for hit in hits:
                    try:
                        aa_seq = hit.getSequence()
                        base_seq = aa_seq.toUnmodifiedString().upper()

                        # Check if base sequence matches any target
                        if base_seq not in target_set:
                            continue

                        # Extract full peptidoform string
                        peptidoform_str = aa_seq.toString()

                        # Get charge
                        charge = hit.getCharge() if hit.getCharge() > 0 else 0

                        # Get score
                        score = hit.getScore() if hit.getScore() else 0.0

                        # Check for id_merge_index on hit level too
                        hit_merge_index = id_merge_index
                        if hit.metaValueExists("id_merge_index"):
                            try:
                                hit_merge_index = int(hit.getMetaValue("id_merge_index"))
                            except (ValueError, TypeError):
                                pass

                        # Re-resolve mzML if hit has different index
                        hit_mzml_filename = mzml_filename
                        hit_metadata_path = mzml_metadata_path
                        if hit_merge_index != id_merge_index and matched_run_id and matched_run_id in run_map:
                            mzml_list = run_map[matched_run_id]
                            full_list = run_map_full.get(matched_run_id, [])
                            if hit_merge_index >= 0 and hit_merge_index < len(mzml_list):
                                hit_mzml_filename = mzml_list[hit_merge_index]
                                if hit_merge_index < len(full_list):
                                    hit_metadata_path = full_list[hit_merge_index]

                        # Resolve to a specific selected mzML path
                        resolved_path = ""
                        if mzml_path_strs and hit_mzml_filename:
                            resolved = resolve_mzml_path(
                                hit_mzml_filename, hit_metadata_path,
                                mzml_path_strs, idxml_path
                            )
                            if resolved:
                                resolved_path = resolved

                        # Create match object
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

                        # Get modification key and create variant key
                        mod_key = get_modification_string(hit)
                        variant_key = f"{base_seq}|{mod_key}" if mod_key else base_seq

                        # Add to variant group
                        if variant_key not in variant_groups:
                            if mod_key:
                                display_label = f"{base_seq} [{mod_key}]"
                            else:
                                display_label = base_seq

                            variant_groups[variant_key] = VariantGroup(
                                base_sequence=base_seq,
                                modification_key=mod_key,
                                display_label=display_label
                            )

                        variant_groups[variant_key].matches.append(match)
                        total_psms_found += 1

                    except Exception as e:
                        logger.debug(f"Error processing peptide hit: {e}")
                        continue

        except Exception as e:
            error_msg = f"Error loading {Path(idxml_path).name}: {e}"
            if error_callback:
                error_callback(error_msg)
            logger.error(error_msg)
            continue

    if progress_callback:
        progress_callback(total_files, total_files, "Search complete")

    # Log summary
    logger.info(f"Found {total_psms_found} PSMs across {len(variant_groups)} variants")

    return variant_groups


def get_unique_peptides_from_idxmls(
    idxml_paths: List[str],
    limit: int = 100,
    progress_callback: Optional[callable] = None
) -> List[str]:
    """
    Extract unique peptide sequences from idXML files.

    Args:
        idxml_paths: List of idXML file paths
        limit: Maximum number of unique peptides to return
        progress_callback: Optional callback for progress

    Returns:
        List of unique base sequences (sorted alphabetically)
    """
    unique_peptides: Set[str] = set()
    total_files = len(idxml_paths)

    for idx, idxml_path in enumerate(idxml_paths):
        if progress_callback:
            progress_callback(idx, total_files, f"Scanning {Path(idxml_path).name}")

        if len(unique_peptides) >= limit:
            break

        try:
            protein_ids = []
            peptide_ids = []
            oms.IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)

            for pep_id in peptide_ids:
                for hit in pep_id.getHits():
                    try:
                        base_seq = hit.getSequence().toUnmodifiedString().upper()
                        unique_peptides.add(base_seq)

                        if len(unique_peptides) >= limit:
                            break
                    except Exception:
                        continue

                if len(unique_peptides) >= limit:
                    break

        except Exception as e:
            logger.warning(f"Error scanning {idxml_path}: {e}")
            continue

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
