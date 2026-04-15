"""Data structures returned by peptide search.

``PeptideMatch`` represents a single PSM with all metadata needed to resolve
its mzML location; ``VariantGroup`` groups matches that share a base sequence
+ modification combination. The pickle helper allows safe serialisation
across Streamlit reloads (where module identity may otherwise differ).
"""
import importlib
from dataclasses import dataclass, field
from typing import List


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
    mod = importlib.import_module("mspecviz.search.models")
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
