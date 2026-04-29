"""
Shared ion-related utilities for MS/MS spectrum handling.
"""
from typing import Dict, Optional

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
