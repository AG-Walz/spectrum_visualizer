"""Peptide search across idXML files.

Re-exports the public API of both submodules so callers can use
``from mspecviz.search import ...`` without caring which file a name lives in.
"""
from mspecviz.search.peptide_search import (
    get_unique_peptides_from_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)
from mspecviz.search.idxml_workers import count_psms_per_peptide

__all__ = [
    "get_unique_peptides_from_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
    "count_psms_per_peptide",
]