"""
Spectrum loading from mzML files.

Provides memory-efficient indexed access to spectra using PyOpenMS.
Uses pre-resolved mzML mappings from the search phase.
"""
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable

import numpy as np
import pyopenms as oms

from core.peptide_search import PeptideMatch, VariantGroup

logger = logging.getLogger(__name__)


@dataclass
class LoadingProgress:
    """Progress tracking object for spectrum loading."""
    files_done: int = 0
    total_files: int = 0
    spectra_loaded: int = 0
    spectra_expected: int = 0
    spectra_not_found: int = 0
    current_file: str = ""
    current_file_spectra: int = 0
    current_file_total: int = 0
    errors: List[str] = field(default_factory=list)

    def add_error(self, error: str):
        """Add an error message."""
        self.errors.append(error)

    def file_complete(self):
        """Mark current file as complete."""
        self.files_done += 1

    def spectrum_loaded(self):
        """Mark a spectrum as successfully loaded."""
        self.spectra_loaded += 1
        self.current_file_spectra += 1

    def spectrum_not_found(self):
        """Mark a spectrum as not found."""
        self.spectra_not_found += 1

    def start_file(self, filename: str, expected_spectra: int):
        """Start processing a new file."""
        self.current_file = filename
        self.current_file_spectra = 0
        self.current_file_total = expected_spectra

    @property
    def total_processed(self) -> int:
        """Total spectra processed (found + not found)."""
        return self.spectra_loaded + self.spectra_not_found

    @property
    def progress_fraction(self) -> float:
        """Overall progress as a fraction 0-1."""
        if self.spectra_expected == 0:
            return 0.0
        return self.total_processed / self.spectra_expected


@dataclass
class LoadedSpectrum:
    """A loaded spectrum with metadata."""
    spectrum_id: str
    mz: np.ndarray
    intensity: np.ndarray
    precursor_mz: float
    precursor_charge: int
    peptidoform: str
    base_sequence: str
    source_file: str
    is_reference: bool = False
    rt: float = 0.0
    match: Optional[PeptideMatch] = None


def build_mzml_to_matches_map(
    variant_groups: Dict[str, VariantGroup],
    available_mzml_names: set
) -> Dict[str, List[Tuple[PeptideMatch, str]]]:
    """
    Build a mapping of mzML filename -> list of (match, variant_key) tuples.

    Uses pre-resolved mzml_filename from PeptideMatch.

    Args:
        variant_groups: Dict of variant groups with matches
        available_mzml_names: Set of available mzML filenames (lowercase)

    Returns:
        Dict mapping mzML filename -> list of (match, variant_key) tuples
    """
    mzml_to_matches: Dict[str, List[Tuple[PeptideMatch, str]]] = {}

    for variant_key, vg in variant_groups.items():
        for match in vg.matches:
            if not match.mzml_filename:
                continue

            # Check if this mzML is in our available files (case-insensitive)
            mzml_lower = match.mzml_filename.lower()
            if mzml_lower in available_mzml_names:
                if mzml_lower not in mzml_to_matches:
                    mzml_to_matches[mzml_lower] = []
                mzml_to_matches[mzml_lower].append((match, variant_key))

    return mzml_to_matches


def open_indexed_mzml(mzml_path: str) -> Tuple[oms.OnDiscMSExperiment, Dict[str, int]]:
    """
    Open an mzML file with indexed access and build ID lookup.

    Args:
        mzml_path: Path to mzML file

    Returns:
        Tuple of (OnDiscMSExperiment, native_id_to_index_dict)
    """
    indexed_exp = oms.OnDiscMSExperiment()
    indexed_exp.openFile(str(mzml_path))

    n_spectra = indexed_exp.getNrSpectra()
    id_to_index = {}

    # Detect the ID format from first spectrum
    id_format = "unknown"
    if n_spectra > 0:
        first_id = indexed_exp.getSpectrum(0).getNativeID()
        if "index=" in first_id:
            id_format = "index"
        elif "scan=" in first_id:
            id_format = "scan"
        elif "spectrum=" in first_id:
            id_format = "spectrum"

    id_to_index["__format__"] = id_format
    id_to_index["__n_spectra__"] = n_spectra

    # Build lookup map for all spectra
    for i in range(n_spectra):
        try:
            spec = indexed_exp.getSpectrum(i)
            native_id = spec.getNativeID()
            id_to_index[native_id] = i

            # Also store by extracted number for flexibility
            num_match = re.search(r'(\d+)', native_id)
            if num_match:
                id_to_index[f"_num_{num_match.group(1)}"] = i
        except Exception:
            continue

    return indexed_exp, id_to_index


def clean_spectrum_id(spectrum_id: str) -> str:
    """Clean a spectrum ID for lookup."""
    if not spectrum_id:
        return ""

    if isinstance(spectrum_id, bytes):
        spectrum_id = spectrum_id.decode()

    clean_id = spectrum_id.strip()

    # If it contains "controllerType" (Thermo format), keep it
    if "controllerType" in clean_id:
        return clean_id

    # For simple scan=N format, ensure clean extraction
    for tag in ["scan=", "index=", "spectrum="]:
        if tag in clean_id:
            idx = clean_id.find(tag)
            rest = clean_id[idx:]
            parts = rest.split()
            return parts[0]

    return clean_id


def find_spectrum_index(
    spectrum_id: str,
    id_to_index: Dict[str, int],
    n_spectra: int
) -> Optional[int]:
    """
    Find spectrum index using the pre-built lookup.

    Args:
        spectrum_id: Native ID to find
        id_to_index: Dictionary mapping native IDs to indices
        n_spectra: Total number of spectra in file

    Returns:
        Spectrum index or None if not found
    """
    clean_id = clean_spectrum_id(spectrum_id)
    if not clean_id:
        return None

    id_format = id_to_index.get("__format__", "unknown")

    # Strategy 1: Direct lookup
    if clean_id in id_to_index:
        return id_to_index[clean_id]

    # Strategy 2: Extract number and try direct index access
    scan_match = re.search(r'(?:scan|index|spectrum)=(\d+)', clean_id)
    if scan_match:
        num = int(scan_match.group(1))

        # If format is "index", the number IS the index
        if id_format == "index" and 0 <= num < n_spectra:
            return num

        # Try numeric lookup
        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

    # Strategy 3: Try just the number
    num_match = re.search(r'(\d+)', clean_id)
    if num_match:
        num = int(num_match.group(1))
        num_key = f"_num_{num}"
        if num_key in id_to_index:
            return id_to_index[num_key]

        if id_format == "index" and 0 <= num < n_spectra:
            return num

    return None


def load_matched_spectra(
    variant_groups: Dict[str, VariantGroup],
    mzml_paths: List[Path],
    reference_substring: str = "",
    progress_callback: Optional[Callable[[LoadingProgress], None]] = None,
    error_callback: Optional[Callable[[str], None]] = None
) -> Dict[str, VariantGroup]:
    """
    Load spectra for all matches in variant groups.

    Uses pre-resolved mzML filenames from the search phase.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup with matches
        mzml_paths: List of available mzML paths
        reference_substring: Substring to identify reference files
        progress_callback: Optional callback(LoadingProgress) for progress updates
        error_callback: Optional callback(error_message) for immediate error reporting

    Returns:
        Updated variant_groups with loaded spectra
    """
    # Build available mzML lookup (case-insensitive)
    mzml_path_map = {p.name.lower(): p for p in mzml_paths}
    available_mzml_names = set(mzml_path_map.keys())

    # Build mapping using pre-resolved filenames
    mzml_to_matches = build_mzml_to_matches_map(variant_groups, available_mzml_names)

    # Calculate totals for progress
    total_expected = sum(len(matches) for matches in mzml_to_matches.values())

    # Also count unresolved matches
    total_unresolved = 0
    for vg in variant_groups.values():
        for match in vg.matches:
            if not match.mzml_filename or match.mzml_filename.lower() not in available_mzml_names:
                total_unresolved += 1

    if total_unresolved > 0:
        error_msg = f"{total_unresolved} PSMs could not be mapped to available mzML files"
        if error_callback:
            error_callback(error_msg)
        logger.warning(error_msg)

    if not mzml_to_matches:
        error_msg = "No mzML files could be mapped to matches"
        if error_callback:
            error_callback(error_msg)
        logger.warning(error_msg)
        return variant_groups

    # Initialize progress
    progress = LoadingProgress(
        total_files=len(mzml_to_matches),
        spectra_expected=total_expected
    )

    # Track loaded spectra by variant
    loaded_spectra: Dict[str, List[LoadedSpectrum]] = {k: [] for k in variant_groups.keys()}

    # Process each mzML file
    for mzml_lower, matches_for_file in mzml_to_matches.items():
        mzml_path = mzml_path_map[mzml_lower]
        is_reference = reference_substring and reference_substring.lower() in mzml_path.name.lower()

        progress.start_file(mzml_path.name, len(matches_for_file))

        if progress_callback:
            progress_callback(progress)

        try:
            # Open indexed file
            indexed_exp, id_to_index = open_indexed_mzml(str(mzml_path))
            n_spectra = id_to_index.get("__n_spectra__", indexed_exp.getNrSpectra())
            logger.info(f"Opened {mzml_path.name} ({n_spectra} spectra, {len(matches_for_file)} to load)")

            # Load each matched spectrum
            for match, variant_key in matches_for_file:
                try:
                    idx = find_spectrum_index(match.spectrum_id, id_to_index, n_spectra)

                    if idx is None:
                        progress.spectrum_not_found()
                        error_msg = f"Spectrum not found: {match.spectrum_id} in {mzml_path.name}"
                        progress.add_error(error_msg)
                        if error_callback:
                            error_callback(error_msg)
                        continue

                    spectrum = indexed_exp.getSpectrum(idx)

                    # Extract data
                    mz_array = np.array([p.getMZ() for p in spectrum])
                    intensity_array = np.array([p.getIntensity() for p in spectrum])

                    # Get precursor info
                    precursor_mz = 0.0
                    precursor_charge = match.charge
                    if spectrum.getPrecursors():
                        prec = spectrum.getPrecursors()[0]
                        precursor_mz = prec.getMZ()
                        if prec.getCharge() > 0:
                            precursor_charge = prec.getCharge()

                    # Create LoadedSpectrum
                    loaded = LoadedSpectrum(
                        spectrum_id=match.spectrum_id,
                        mz=mz_array,
                        intensity=intensity_array,
                        precursor_mz=precursor_mz,
                        precursor_charge=precursor_charge,
                        peptidoform=match.peptidoform,
                        base_sequence=match.base_sequence,
                        source_file=mzml_path.name,
                        is_reference=is_reference,
                        rt=spectrum.getRT() if spectrum.getRT() else 0.0,
                        match=match
                    )

                    loaded_spectra[variant_key].append(loaded)
                    progress.spectrum_loaded()

                except Exception as e:
                    progress.spectrum_not_found()
                    error_msg = f"Error loading {match.spectrum_id}: {e}"
                    progress.add_error(error_msg)
                    if error_callback:
                        error_callback(error_msg)
                    continue

            progress.file_complete()

            # Update progress after each file
            if progress_callback:
                progress_callback(progress)

        except Exception as e:
            error_msg = f"Error processing {mzml_path.name}: {e}"
            progress.add_error(error_msg)
            if error_callback:
                error_callback(error_msg)
            logger.error(error_msg)
            progress.file_complete()
            if progress_callback:
                progress_callback(progress)
            continue

    # Assign loaded spectra to variant groups
    for variant_key, spectra in loaded_spectra.items():
        if variant_key in variant_groups:
            vg = variant_groups[variant_key]
            for spec in spectra:
                if spec.is_reference:
                    vg.reference_spectra.append(spec)
                else:
                    vg.query_spectra.append(spec)

    # Log summary
    logger.info(f"Loaded {progress.spectra_loaded} spectra, {progress.spectra_not_found} not found, across {progress.files_done} files")

    return variant_groups
