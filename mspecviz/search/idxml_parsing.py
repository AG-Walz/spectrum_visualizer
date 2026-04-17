"""idXML metadata + peptidoform parsing helpers.

Low-level pyopenms / psm_utils facade: spectra-data extraction and compact
modification-key formatting. The search workers and orchestration layers in
this package import from here.
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
            raw_value = protein_id.getMetaValue("spectra_data")
            raw_list = raw_value if isinstance(raw_value, list) else [raw_value]
            for entry in raw_list:
                full_path = decode_metadata_value(entry)
                mzml_full_paths.append(full_path)
                mzml_files.append(Path(full_path).name)

        # Fallback metadata keys (single-file idXMLs, older exporters)
        if not mzml_files:
            for key in ("CometAdapter:1:in", "mzML_file", "spectrum_reference"):
                if protein_id.metaValueExists(key):
                    decoded = decode_metadata_value(protein_id.getMetaValue(key))
                    if decoded:
                        mzml_full_paths = [decoded]
                        mzml_files = [Path(decoded).name]
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
        for mod in peptidoform_obj.properties["n_term"]:
            mods.append(f"NTerm[{abbreviate_mod(mod.name)}]")
    for position, (amino_acid, residue_mods) in enumerate(peptidoform_obj.parsed_sequence):
        if residue_mods:
            for mod in residue_mods:
                mods.append(f"{amino_acid}{position + 1}[{abbreviate_mod(mod.name)}]")
    if peptidoform_obj.properties.get("c_term"):
        for mod in peptidoform_obj.properties["c_term"]:
            mods.append(f"CTerm[{abbreviate_mod(mod.name)}]")
    return "|".join(mods) if mods else ""
