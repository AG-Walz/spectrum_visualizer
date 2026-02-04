"""
Mirror plot generation for MS/MS spectrum comparison.

Creates publication-ready mirror plots using Matplotlib.
"""
import logging
from io import BytesIO
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default ion colors
ION_COLORS = {
    "a": "#FFCC00",
    "b": "#0066CC",
    "c": "#00CC66",
    "x": "#CC6600",
    "y": "#CC0000",
    "z": "#6600CC",
    "precursor": "#8B008B",
    "immonium": "#00CED1",
    "unannotated": "#AAAAAA"
}


def get_ion_color(annotation: str, color_map: Dict[str, str] = None) -> str:
    """
    Get color for an ion annotation.

    Args:
        annotation: Ion annotation string
        color_map: Optional custom color mapping

    Returns:
        Hex color string
    """
    if color_map is None:
        color_map = ION_COLORS

    if not annotation or annotation == "unannotated":
        return color_map.get("unannotated", "#AAAAAA")

    annotation_lower = annotation.lower()

    # Check for precursor
    if annotation_lower.startswith("[m") or "precursor" in annotation_lower:
        return color_map.get("precursor", "#8B008B")

    # Check for standard ions
    for ion in ["a", "b", "c", "x", "y", "z"]:
        if annotation_lower.startswith(ion):
            return color_map.get(ion, "#666666")

    # Check for immonium
    if "immonium" in annotation_lower:
        return color_map.get("immonium", "#00CED1")

    return color_map.get("unannotated", "#AAAAAA")


def create_mirror_plot(
    upper_spectrum: pd.DataFrame,
    lower_spectrum: pd.DataFrame,
    title: str = "Mirror Plot",
    upper_label: str = "Query",
    lower_label: str = "Reference",
    color_map: Dict[str, str] = None,
    figsize: Tuple[int, int] = (12, 8),
    dpi: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7
):
    """
    Create a publication-ready mirror plot.

    Args:
        upper_spectrum: DataFrame with columns [mz, intensity, annotation]
        lower_spectrum: DataFrame with columns [mz, intensity, annotation]
        title: Plot title
        upper_label: Label for upper spectrum (query)
        lower_label: Label for lower spectrum (reference)
        color_map: Custom ion color mapping
        figsize: Figure size (width, height)
        dpi: Resolution for export
        show_annotations: Whether to show ion labels
        annotation_fontsize: Font size for annotations

    Returns:
        matplotlib.figure.Figure object
    """
    # Import matplotlib with non-interactive backend
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure

    if color_map is None:
        color_map = ION_COLORS

    # Create figure
    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    # Process upper spectrum (positive intensities)
    if not upper_spectrum.empty and len(upper_spectrum) > 0:
        upper_df = upper_spectrum.copy()

        # Ensure annotation column exists
        if 'annotation' not in upper_df.columns:
            upper_df['annotation'] = 'unannotated'

        # Normalize
        max_int = upper_df["intensity"].max()
        if max_int > 0:
            upper_df["intensity"] = upper_df["intensity"] / max_int

        # Plot each peak
        for _, peak in upper_df.iterrows():
            annotation = peak.get('annotation', 'unannotated')
            color = get_ion_color(annotation, color_map)
            ax.vlines(peak['mz'], 0, peak['intensity'],
                     colors=color, linewidth=1.5, alpha=0.8)

            # Add annotation label
            if show_annotations and annotation and annotation != 'unannotated':
                ax.text(peak['mz'], peak['intensity'] + 0.03, annotation,
                       fontsize=annotation_fontsize, color=color,
                       ha='center', va='bottom', rotation=0)

    # Process lower spectrum (negative/mirrored intensities)
    if not lower_spectrum.empty and len(lower_spectrum) > 0:
        lower_df = lower_spectrum.copy()

        # Ensure annotation column exists
        if 'annotation' not in lower_df.columns:
            lower_df['annotation'] = 'unannotated'

        # Normalize and mirror
        max_int = lower_df["intensity"].max()
        if max_int > 0:
            lower_df["intensity"] = (lower_df["intensity"] / max_int) * -1

        # Plot each peak
        for _, peak in lower_df.iterrows():
            annotation = peak.get('annotation', 'unannotated')
            color = get_ion_color(annotation, color_map)
            ax.vlines(peak['mz'], 0, peak['intensity'],
                     colors=color, linewidth=1.5, alpha=0.8)

            # Add annotation label
            if show_annotations and annotation and annotation != 'unannotated':
                ax.text(peak['mz'], peak['intensity'] - 0.03, annotation,
                       fontsize=annotation_fontsize, color=color,
                       ha='center', va='top', rotation=0)

    # Set axis limits
    all_mz = []
    if not upper_spectrum.empty:
        all_mz.extend(upper_spectrum["mz"].tolist())
    if not lower_spectrum.empty:
        all_mz.extend(lower_spectrum["mz"].tolist())

    if all_mz:
        mz_min = min(all_mz) - 50
        mz_max = max(all_mz) + 50
        ax.set_xlim(mz_min, mz_max)

    ax.set_ylim(-1.15, 1.15)

    # Add center line
    ax.axhline(y=0, color='black', linewidth=1, alpha=0.5)

    # Add labels at the sides
    ax.text(0.02, 0.95, upper_label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left',
            color='#333333')
    ax.text(0.02, 0.05, lower_label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='left',
            color='#333333')

    # Styling
    ax.set_xlabel("m/z", fontsize=12, fontweight='bold')
    ax.set_ylabel("Relative Intensity", fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Custom y-axis labels
    ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
    ax.set_yticklabels(['1.0', '0.5', '0', '0.5', '1.0'])

    fig.tight_layout()

    return fig


def create_legend_figure(
    color_map: Dict[str, str] = None,
    ion_types: List[str] = None,
    figsize: Tuple[float, float] = (4, 2),
    dpi: int = 150
):
    """
    Create a standalone legend figure for ion colors.

    Args:
        color_map: Ion color mapping
        ion_types: Ion types to include (default: all)
        figsize: Figure size
        dpi: Resolution

    Returns:
        matplotlib Figure with legend
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    if color_map is None:
        color_map = ION_COLORS

    if ion_types is None:
        ion_types = list(color_map.keys())

    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    ax.axis('off')

    # Create legend patches
    patches = []
    for ion_type in ion_types:
        if ion_type in color_map:
            label = ion_type if ion_type != "unannotated" else "unassigned"
            patches.append(Patch(facecolor=color_map[ion_type], label=label))

    ax.legend(handles=patches, loc='center', ncol=min(4, len(patches)),
              frameon=False, fontsize=10)

    fig.tight_layout()
    return fig


def figure_to_bytes(fig, format: str = 'png') -> bytes:
    """
    Convert a matplotlib figure to bytes.

    Args:
        fig: matplotlib Figure
        format: Output format ('png' or 'svg')

    Returns:
        Bytes of the image
    """
    buf = BytesIO()
    fig.savefig(buf, format=format, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    return buf.read()


def figure_to_base64(fig, format: str = 'png') -> str:
    """
    Convert a matplotlib figure to base64 encoded string.

    Useful for embedding in HTML.
    """
    import base64
    img_bytes = figure_to_bytes(fig, format)
    return base64.b64encode(img_bytes).decode('utf-8')
