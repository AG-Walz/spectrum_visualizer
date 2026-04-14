"""
idXML metadata + peptidoform parsing helpers.

Extracted from ``peptide_search`` so the low-level pyopenms/psm_utils facade
lives in a focused module. ``peptide_search`` keeps the public search API
(``search_peptides_in_idxmls`` etc.), variant grouping, and the worker
function; it imports ``extract_spectra_data_map`` and
``_mod_key_from_peptidoform`` from here.
"""
from pathlib import Path
from typing import Dict, List, Tuple

# Compact abbreviations for common post-translational modifications
MOD_ABBREVIATIONS = {
    "Oxidation": "Ox",
    "Phospho": "Ph",
    "Carbamidomethyl": "CAM",
    "Acetyl": "Ac",
    "Deamidated": "Deam",
    "Methylation": "Me",
    "Dimethyl": "DiMe",
    "Ubiquitin": "GG",
    "Pyro-carbamidomethyl": "Pyro-CAM",
    "Pyro-glu": "Pyro",
}


def abbreviate_mod(name: str) -> str:
    """Return short abbreviation for a modification, or the full name if unknown."""
    return MOD_ABBREVIATIONS.get(name, name)


def decode_metadata_value(value) -> str:
    """Decode metadata value that might be bytes, list of bytes, or string."""
    if isinstance(value, list) and len(value) > 0:
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value) if value else ""


def extract_spectra_data_map(protein_ids: List) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Extract spectra_data mapping from ProteinIdentification objects.

    Returns:
        Tuple of:
        - run_map: identifier -> list of mzML bare filenames
        - run_map_full: identifier -> list of full metadata paths (as stored in idXML)
    """
    run_map: Dict[str, List[str]] = {}
    run_map_full: Dict[str, List[str]] = {}

    for protein_id in protein_ids:
        identifier = protein_id.getIdentifier()
        mzml_files: List[str] = []
        mzml_full_paths: List[str] = []

        # Standard path: spectra_data (merged idXML)
        if protein_id.metaValueExists("spectra_data"):
            raw_val = protein_id.getMetaValue("spectra_data")
            raw_list = raw_val if isinstance(raw_val, list) else [raw_val]
            for v in raw_list:
                full = decode_metadata_value(v)
                mzml_full_paths.append(full)
                mzml_files.append(Path(full).name)

        # Fallback metadata keys (single-file idXMLs, older exporters)
        if not mzml_files:
            for key in ("CometAdapter:1:in", "mzML_file", "spectrum_reference"):
                if protein_id.metaValueExists(key):
                    val = decode_metadata_value(protein_id.getMetaValue(key))
                    if val:
                        mzml_full_paths = [val]
                        mzml_files = [Path(val).name]
                        break

        if mzml_files:
            run_map[identifier] = mzml_files
            run_map_full[identifier] = mzml_full_paths

    return run_map, run_map_full


def mod_key_from_peptidoform(peptidoform_obj) -> str:
    """
    Build compact modification key from a psm_utils Peptidoform object.

    Returns format like "M5[Ox]" or "NTerm[Ac]|M5[Ox]|S3[Ph]" for multiple mods.
    """
    mods = []
    if peptidoform_obj.properties.get("n_term"):
        for m in peptidoform_obj.properties["n_term"]:
            mods.append(f"NTerm[{abbreviate_mod(m.name)}]")
    for i, (aa, residue_mods) in enumerate(peptidoform_obj.parsed_sequence):
        if residue_mods:
            for m in residue_mods:
                mods.append(f"{aa}{i + 1}[{abbreviate_mod(m.name)}]")
    if peptidoform_obj.properties.get("c_term"):
        for m in peptidoform_obj.properties["c_term"]:
            mods.append(f"CTerm[{abbreviate_mod(m.name)}]")
    return "|".join(mods) if mods else ""
