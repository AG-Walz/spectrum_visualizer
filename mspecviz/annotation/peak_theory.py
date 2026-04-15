"""Theoretical fragment generation for peptidoforms.

Converts ProForma peptidoforms to PyOpenMS ``AASequence`` objects and runs
``TheoreticalSpectrumGenerator`` to produce the candidate fragment table that
:mod:`mspecviz.annotation.peak_matching` then aligns to experimental peaks.
"""
import functools
import logging
from typing import List

import pandas as pd
import pyopenms as oms
from psm_utils import Peptidoform

logger = logging.getLogger(__name__)

DEFAULT_ION_TYPES = ["b", "y"]
ALL_ION_TYPES = ["a", "b", "c", "x", "y", "z", "precursor", "immonium"]


def peptidoform_to_aasequence(peptidoform_str: str) -> oms.AASequence:
    """
    Convert a ProForma peptidoform string to PyOpenMS AASequence.

    Parses via psm_utils Peptidoform and reconstructs in OpenMS notation,
    correctly handling terminal modifications (e.g. [Acetyl]-PEPTIDE becomes
    .(Acetyl)PEPTIDE rather than the invalid (Acetyl)-PEPTIDE).
    """
    clean_peptidoform_str = peptidoform_str.split('/')[0]
    peptidoform_obj = Peptidoform(clean_peptidoform_str)

    # Rebuild in OpenMS parenthesis notation from the parsed representation
    parts = []
    if peptidoform_obj.properties.get("n_term"):
        for m in peptidoform_obj.properties["n_term"]:
            parts.append(f".({m.name})")
    for aa, residue_mods in peptidoform_obj.parsed_sequence:
        parts.append(aa)
        if residue_mods:
            for m in residue_mods:
                parts.append(f"({m.name})")
    if peptidoform_obj.properties.get("c_term"):
        for m in peptidoform_obj.properties["c_term"]:
            parts.append(f".({m.name})")
    oms_str = "".join(parts)

    try:
        return oms.AASequence.fromString(oms_str)
    except Exception as e:
        logger.warning(f"Could not parse peptidoform '{peptidoform_str}': {e}")
        try:
            return oms.AASequence.fromString(peptidoform_obj.sequence)
        except Exception:
            raise ValueError(f"Cannot convert peptidoform: {peptidoform_str}")


@functools.lru_cache(maxsize=1024)
def _generate_fragments_impl(
    peptidoform: str,
    charge: int,
    ion_types: tuple,          # tuple (hashable) — converted from list by caller
    max_fragment_charge: int,
    add_losses: bool,
    add_first_prefix_ion: bool,
) -> pd.DataFrame:
    """
    Cached implementation of theoretical fragment generation.

    Returns the same DataFrame object on repeated calls with identical
    arguments — callers must not mutate the returned DataFrame.
    """
    try:
        aa_seq = peptidoform_to_aasequence(peptidoform)
    except ValueError as e:
        logger.error(f"Failed to generate fragments: {e}")
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Configure theoretical spectrum generator
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()

    ion_types_list = list(ion_types)
    for ion in ["a", "b", "c", "x", "y", "z"]:
        params.setValue(f"add_{ion}_ions", "true" if ion in ion_types_list else "false")

    params.setValue("add_metainfo", "true")
    params.setValue("add_losses", "true" if add_losses else "false")
    params.setValue("add_precursor_peaks", "true")
    params.setValue("add_all_precursor_charges", "true")
    params.setValue("add_first_prefix_ion", "true" if add_first_prefix_ion else "false")
    params.setValue("add_abundant_immonium_ions", "true" if "immonium" in ion_types_list else "false")

    tsg.setParameters(params)

    # Generate spectrum
    theo_spec = oms.MSSpectrum()
    tsg.getSpectrum(theo_spec, aa_seq, 1, max_fragment_charge)

    if theo_spec.size() > 1:
        theo_spec.sortByPosition()

    if theo_spec.size() == 0:
        return pd.DataFrame(columns=["mz", "intensity", "annotation", "charge"])

    # Use batch API — much faster than iterating Peak1D objects
    mzs, ints = theo_spec.get_peaks()

    # Extract annotations from StringDataArrays
    string_arrays = theo_spec.getStringDataArrays()
    int_arrays = theo_spec.getIntegerDataArrays()

    if string_arrays:
        annotations = [n.decode() if isinstance(n, bytes) else str(n) for n in string_arrays[0]]
    else:
        annotations = [""] * theo_spec.size()

    if int_arrays:
        charges = list(int_arrays[0])
    else:
        charges = [1] * theo_spec.size()

    return pd.DataFrame({
        "mz": mzs,
        "intensity": ints,
        "annotation": annotations,
        "charge": charges,
    })


def generate_theoretical_fragments(
    peptidoform: str,
    charge: int = 2,
    ion_types: List[str] = None,
    max_fragment_charge: int = 2,
    add_losses: bool = False,
    add_first_prefix_ion: bool = True,
) -> pd.DataFrame:
    """
    Generate theoretical fragment ions for a peptide.

    Results are cached by (peptidoform, charge, ion_types, add_losses,
    add_first_prefix_ion) using functools.lru_cache — identical calls within
    the same process (e.g. many query spectra for the same peptidoform) pay
    the PyOpenMS cost only once.

    Args:
        peptidoform: Peptidoform string (with or without modifications)
        charge: Precursor charge state
        ion_types: List of ion types to generate (default: ["b", "y"])
        max_fragment_charge: Maximum charge for fragment ions
        add_losses: Whether to add neutral loss peaks (H2O, NH3)
        add_first_prefix_ion: Whether to add first prefix ion (a1, b1)

    Returns:
        DataFrame with columns: mz, intensity (relative), annotation, charge
        Do NOT modify the returned DataFrame — it is shared across callers.
    """
    if ion_types is None:
        ion_types = DEFAULT_ION_TYPES

    return _generate_fragments_impl(
        peptidoform,
        charge,
        tuple(ion_types),
        max_fragment_charge,
        add_losses,
        add_first_prefix_ion,
    )
