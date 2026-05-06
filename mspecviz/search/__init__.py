"""Peptide search across idXML files.

Re-exports the public API of the four submodules so callers can keep using
``from mspecviz.search import ...`` without caring which file a name lives in.
:class:`psm_utils.PSM` is also re-exported under the name ``PSM`` (it replaces
the former local ``PeptideMatch`` dataclass).
"""
from psm_utils import PSM

from mspecviz.search.idxml_parsing import (
    MOD_ABBREVIATIONS,
    abbreviate_modification,
)
from mspecviz.search.models import (
    VariantGroup,
    display_label,
)
from mspecviz.search.orchestration import (
    get_unique_peptides_from_idxmls,
    search_peptides_in_idxmls,
    summarize_search_results,
)
from mspecviz.search.workers import count_psms_per_peptide

__all__ = [
    # psm_utils re-export (replacement for the deleted PeptideMatch)
    "PSM",
    # idxml_parsing
    "MOD_ABBREVIATIONS",
    "abbreviate_modification",
    # models
    "VariantGroup",
    "display_label",
    # orchestration
    "get_unique_peptides_from_idxmls",
    "search_peptides_in_idxmls",
    "summarize_search_results",
    # workers
    "count_psms_per_peptide",
]
