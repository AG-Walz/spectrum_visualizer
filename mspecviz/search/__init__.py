"""Peptide search across idXML files.

Re-exports the public API of the four submodules so callers can keep using
``from mspecviz.search import ...`` without caring which file a name lives in.
"""
from mspecviz.search.idxml_parsing import (
    MOD_ABBREVIATIONS,
    abbreviate_mod,
    decode_metadata_value,
    extract_spectra_data_map,
    mod_key_from_peptidoform,
)
from mspecviz.search.models import (
    PeptideMatch,
    VariantGroup,
    format_inline_sequence,
)
from mspecviz.search.orchestration import (
    get_unique_peptides_from_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)
from mspecviz.search.workers import count_psms_per_peptide

__all__ = [
    # idxml_parsing
    "MOD_ABBREVIATIONS",
    "abbreviate_mod",
    "decode_metadata_value",
    "extract_spectra_data_map",
    "mod_key_from_peptidoform",
    # models
    "PeptideMatch",
    "VariantGroup",
    "format_inline_sequence",
    # workers
    "count_psms_per_peptide",
    # orchestration
    "get_unique_peptides_from_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
]
