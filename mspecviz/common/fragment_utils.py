"""Fragment-ion + peptidoform primitives shared across annotation, plotting, and prediction.

Two concerns live here:
  * Fragment-ion typing/colours — ``get_ion_type``, ``get_ion_color``, ``ION_COLORS``.
  * ProForma ↔ PyOpenMS conversion — ``peptidoform_to_aasequence``.

Lives in ``common/`` so neither ``annotation/`` nor ``spectra/`` has to depend
on the other.
"""
import logging
from typing import Dict, Optional

import pyopenms as oms
from psm_utils import Peptidoform

logger = logging.getLogger(__name__)


ION_COLORS = {
    "a": "#FFCC00",
    "b": "#0066CC",
    "c": "#00CC66",
    "x": "#CC6600",
    "y": "#CC0000",
    "z": "#6600CC",
    "precursor": "#8B008B",
    "immonium": "#00CED1",
    "unannotated": "#AAAAAA",
}


def get_ion_type(annotation: str) -> str:
    """
    Extract ion type from annotation string.
        "b3++" -> "b"
        "y7" -> "y"
        "[M+2H]2+" -> "precursor"
        "unannotated" -> "unannotated"
    """
    if not annotation or annotation == "unannotated":
        return "unannotated"

    annotation_lower = annotation.lower()

    if annotation_lower.startswith("[m") or "precursor" in annotation_lower:
        return "precursor"

    for ion in ["a", "b", "c", "x", "y", "z"]:
        if annotation_lower.startswith(ion) and len(annotation) > 1:
            if annotation[1].isdigit():
                return ion

    if "immonium" in annotation_lower:
        return "immonium"

    return "unannotated"


def get_ion_color(annotation: str, color_map: Optional[Dict[str, str]] = None) -> str:
    """Get color for an ion annotation."""
    if color_map is None:
        color_map = ION_COLORS
    fallback = color_map.get("unannotated", "#AAAAAA")
    return color_map.get(get_ion_type(annotation), fallback)


def peptidoform_to_aasequence(peptidoform_str: str) -> oms.AASequence:
    """Parse a ProForma peptidoform via psm_utils, then rebuild it in PyOpenMS notation.

    Handles terminal modifications correctly (``[Acetyl]-PEPTIDE`` becomes
    ``.(Acetyl)PEPTIDE`` rather than the invalid ``(Acetyl)-PEPTIDE``).
    """
    clean_peptidoform_str = peptidoform_str.split("/")[0]
    peptidoform_obj = Peptidoform(clean_peptidoform_str)

    parts = []
    if peptidoform_obj.properties.get("n_term"):
        for mod in peptidoform_obj.properties["n_term"]:
            parts.append(f".({mod.name})")
    for amino_acid, residue_mods in peptidoform_obj.parsed_sequence:
        parts.append(amino_acid)
        if residue_mods:
            for mod in residue_mods:
                parts.append(f"({mod.name})")
    if peptidoform_obj.properties.get("c_term"):
        for mod in peptidoform_obj.properties["c_term"]:
            parts.append(f".({mod.name})")
    oms_str = "".join(parts)

    try:
        return oms.AASequence.fromString(oms_str)
    except RuntimeError as e:
        logger.warning(f"Could not parse peptidoform '{peptidoform_str}': {e}")
    try:
        return oms.AASequence.fromString(peptidoform_obj.sequence)
    except RuntimeError:
        raise ValueError(f"Cannot convert peptidoform: {peptidoform_str}")
