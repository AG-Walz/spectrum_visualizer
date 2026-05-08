"""Peptide search across idXML files.

Re-exports the public API of both submodules so callers can keep using
``from mspecviz.search import ...`` without calling the exacth files.
"""
from psm_utils import PSM, PSMList

from mspecviz.search.peptide_search import (
    get_unique_peptides_from_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)
from mspecviz.search.idxml_workers import count_psms_per_peptide

__all__ = [
    # psm_utils re-exports (replacements for the deleted PeptideMatch / VariantGroup)
    "PSM",
    "PSMList",
    # orchestration
    "get_unique_peptides_from_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
    # workers
    "count_psms_per_peptide",
]
