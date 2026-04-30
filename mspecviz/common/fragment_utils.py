"""Fragment-ion typing + colour palette shared by annotation and plotting.

In ``common/`` so neither ``annotation/`` nor ``plotting/`` has to depend on
the other.
"""
# Per spectrum_utils' SUPPORTED_IONS contract ("?abcxyzIm_prf"), the leading
# character of FragmentAnnotation.ion_type identifies the ion class. Koina
# annotation strings follow the same convention ("y3+1", "b5++"). We keep the
# six fragment letters as their letter so user selections like ['b', 'y']
# match the dataframe; uppercase 'I' and lowercase 'p' map to readable buckets.
_BUCKET_FROM_PREFIX = {
    "a": "a", "b": "b", "c": "c", "x": "x", "y": "y", "z": "z",
    "I": "immonium", "p": "precursor",
}

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
    """Map a spectrum_utils ion_type or Koina annotation to its ``ION_COLORS`` bucket."""
    return _BUCKET_FROM_PREFIX.get(annotation[:1], "unannotated") if annotation else "unannotated"


def get_ion_color(annotation: str) -> str:
    """Return the colour for an ion annotation; falls back to the 'unannotated' colour."""
    return ION_COLORS[get_ion_type(annotation)]
