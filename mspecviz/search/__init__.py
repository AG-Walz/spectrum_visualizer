"""Peptide search across idXML files.

Re-exports the public API of both submodules so callers can use
``from mspecviz.search import ...`` without caring which file a name lives in.
"""
from mspecviz.search.peptide_search import (
    count_psms_in_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)

__all__ = [
    "count_psms_in_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
]