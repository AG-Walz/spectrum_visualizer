"""Modification name abbreviations for mspecviz peptide rendering."""

# Compact abbreviations for common post-translational modifications. Used by
# :func:`mspecviz.search.models.display_label` to abbreviate ProForma names
# (``Oxidation`` → ``Ox``) for inline peptide rendering in tables and plots.
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


def abbreviate_modification(name: str) -> str:
    """Return short abbreviation for a modification, or the full name if unknown."""
    return MOD_ABBREVIATIONS.get(name, name)
