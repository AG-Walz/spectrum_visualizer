"""
Shared ion-related utilities for MS/MS spectrum handling.

Used by src.peak_annotation and src.mirror_plot.
"""
from typing import Dict

# Default ion colors (same as MSpecViz)
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

    Examples:
        "b3++" -> "b"
        "y7" -> "y"
        "[M+2H]2+" -> "precursor"
        "unannotated" -> "unannotated"
    """
    if not annotation or annotation == "unannotated":
        return "unannotated"

    annotation_lower = annotation.lower()

    # Check for precursor
    if annotation_lower.startswith("[m") or "precursor" in annotation_lower:
        return "precursor"

    # Check for standard ions
    for ion in ["a", "b", "c", "x", "y", "z"]:
        if annotation_lower.startswith(ion) and len(annotation) > 1:
            next_char = annotation[1]
            if next_char.isdigit():
                return ion

    # Check for immonium
    if "immonium" in annotation_lower or annotation.startswith("I"):
        return "immonium"

    return "unannotated"


def get_ion_color(annotation: str, color_map: Dict[str, str] = None) -> str:
    """Get color for an ion annotation."""
    if color_map is None:
        color_map = ION_COLORS

    if not annotation or annotation == "unannotated":
        return color_map.get("unannotated", "#AAAAAA")

    annotation_lower = annotation.lower()

    if annotation_lower.startswith("[m") or "precursor" in annotation_lower:
        return color_map.get("precursor", "#8B008B")

    for ion in ["a", "b", "c", "x", "y", "z"]:
        if annotation_lower.startswith(ion):
            return color_map.get(ion, "#666666")

    if "immonium" in annotation_lower:
        return color_map.get("immonium", "#00CED1")

    return color_map.get("unannotated", "#AAAAAA")
