"""Peptide-search public API (consumed by ``cli.main`` and ``streamlit_app.views.search_load``)."""
from core.search.peptide_search import (
    count_psms_in_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)

__all__ = [
    "count_psms_in_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
]