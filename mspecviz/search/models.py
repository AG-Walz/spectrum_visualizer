"""Peptide-search data structures (``PeptideMatch``, ``VariantGroup``) + pickle-safe helper."""
import importlib
from dataclasses import dataclass, field
from typing import List


def _parse_modification_entry(entry: str):
    """Parse one ``M5[Ox]`` / ``NTerm[Ac]`` token into ``(prefix, abbrev)`` or ``None``."""
    entry = entry.strip()
    if not entry:
        return None
    bracket_start = entry.find("[")
    bracket_end = entry.find("]")
    if bracket_start < 0 or bracket_end < 0:
        return None
    return entry[:bracket_start], entry[bracket_start:bracket_end + 1]


def _parse_modification_key(modification_key: str):
    """Split ``modification_key`` into ``(prefix_mods, suffix_mods, position→[abbrev])``."""
    insertions: dict[int, list[str]] = {}
    prefix_mods: list[str] = []
    suffix_mods: list[str] = []

    for entry in modification_key.split("|"):
        parsed = _parse_modification_entry(entry)
        if parsed is None:
            continue
        prefix, abbrev = parsed

        if prefix == "NTerm":
            prefix_mods.append(abbrev)
        elif prefix == "CTerm":
            suffix_mods.append(abbrev)
        else:
            pos_str = "".join(ch for ch in prefix if ch.isdigit())
            if pos_str:
                insertions.setdefault(int(pos_str), []).append(abbrev)

    return prefix_mods, suffix_mods, insertions


def format_inline_sequence(base_sequence: str, modification_key: str) -> str:
    """Splice ``modification_key`` (``"M5[Ox]|NTerm[Ac]"``) into ``base_sequence``."""
    if not modification_key:
        return base_sequence

    prefix_mods, suffix_mods, insertions = _parse_modification_key(modification_key)

    parts: list[str] = list(prefix_mods)
    for i, aa in enumerate(base_sequence):
        parts.append(aa)
        if (pos := i + 1) in insertions:
            parts.extend(insertions[pos])
    parts.extend(suffix_mods)
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
