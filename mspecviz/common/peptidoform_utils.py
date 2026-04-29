"""ProForma ↔ PyOpenMS peptidoform conversion.

Used by Koina prediction (UniMod string for the model input) and by anything
else that needs to hand a parsed peptidoform to PyOpenMS. Lives in ``common/``
so neither ``annotation/`` nor ``spectra/`` has to depend on the other.
"""
import logging

import pyopenms as oms
from psm_utils import Peptidoform

logger = logging.getLogger(__name__)


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
