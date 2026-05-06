"""Peptide-search data structures.

PSM data uses :class:`psm_utils.PSM` directly. mspecviz-specific resolution
metadata (``idxml_file``, ``mzml_resolved_path``, ``idxml_resolved_path``)
is stuffed into ``PSM.metadata`` by ``mspecviz.search.workers``.

This module owns the variant-group container and the display-label renderer
that uses :data:`MOD_ABBREVIATIONS` to convert long modification names
(``Oxidation``) to mspecviz's short form (``Ox``).
"""
from dataclasses import dataclass, field
from psm_utils import Peptidoform
from mspecviz.search.idxml_parsing import MOD_ABBREVIATIONS


def display_label(peptidoform: Peptidoform) -> str:
    """Render a peptidoform as an inline display label using mspecviz mod abbreviations.

    e.g. ``Peptidoform("PEPTM[Oxidation]IDEK/2")`` → ``"PEPTM[Ox]IDEK"``.
    Re-parses ``peptidoform.proforma`` so the caller's object isn't mutated.
    """
    abbreviated = Peptidoform(peptidoform.proforma)
    abbreviated.rename_modifications(MOD_ABBREVIATIONS)
    return abbreviated.modified_sequence


@dataclass
class VariantGroup:
    """Group of PSMs for one peptide variant (sequence + modification state)."""
    base_sequence: str
    modification_key: str     # canonical modified_sequence (or "" for unmodified)
    display_label: str        # human-readable label with abbreviated mod names
    matches: list = field(default_factory=list)  # list[psm_utils.PSM]
