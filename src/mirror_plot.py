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
    annotation_fontsize: int = 7,
    metrics_text: Optional[str] = None,
    normalize: bool = True,
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
        metrics_text: Optional multi-line string of comparison metrics to
            render in the top-right corner of the plot

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

    # Compute annotation text offset (3% of y-range)
    if normalize:
        ann_offset = 0.03
    else:
        upper_max = upper_spectrum["intensity"].max() if not upper_spectrum.empty else 1
        lower_max = lower_spectrum["intensity"].max() if not lower_spectrum.empty else 1
        ann_offset = max(upper_max, lower_max) * 0.03

    # Process upper spectrum (positive intensities)
    if not upper_spectrum.empty and len(upper_spectrum) > 0:
        upper_df = upper_spectrum.copy()

        # Ensure annotation column exists
        if 'annotation' not in upper_df.columns:
            upper_df['annotation'] = 'unannotated'

        # Normalize if requested
        if normalize:
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
                ax.text(peak['mz'], peak['intensity'] + ann_offset, annotation,
                       fontsize=annotation_fontsize, color=color,
                       ha='center', va='bottom', rotation=0)

    # Process lower spectrum (negative/mirrored intensities)
    if not lower_spectrum.empty and len(lower_spectrum) > 0:
        lower_df = lower_spectrum.copy()

        # Ensure annotation column exists
        if 'annotation' not in lower_df.columns:
            lower_df['annotation'] = 'unannotated'

        # Normalize and mirror (or just mirror for raw mode)
        if normalize:
            max_int = lower_df["intensity"].max()
            if max_int > 0:
                lower_df["intensity"] = (lower_df["intensity"] / max_int) * -1
        else:
            lower_df["intensity"] = lower_df["intensity"] * -1

        # Plot each peak
        for _, peak in lower_df.iterrows():
            annotation = peak.get('annotation', 'unannotated')
            color = get_ion_color(annotation, color_map)
            ax.vlines(peak['mz'], 0, peak['intensity'],
                     colors=color, linewidth=1.5, alpha=0.8)

            # Add annotation label
            if show_annotations and annotation and annotation != 'unannotated':
                ax.text(peak['mz'], peak['intensity'] - ann_offset, annotation,
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

    if normalize:
        ax.set_ylim(-1.2, 1.2)
    else:
        # Auto-scale for raw intensities
        upper_max = upper_spectrum["intensity"].max() if not upper_spectrum.empty else 1
        lower_max = lower_spectrum["intensity"].max() if not lower_spectrum.empty else 1
        y_max = max(upper_max, lower_max) * 1.2
        ax.set_ylim(-y_max, y_max)

    # Add center line
    ax.axhline(y=0, color='black', linewidth=1, alpha=0.5)

    # Add labels at the sides - positioned above/below intensity 1.0 to avoid overlap
    # Upper label at top of plot area (above intensity 1.0)
    ax.text(0.02, 0.99, upper_label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left',
            color='#333333',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.8))
    # Lower label at bottom of plot area (below intensity -1.0)
    ax.text(0.02, 0.01, lower_label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='bottom', ha='left',
            color='#333333',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.8))

    # Metrics overlay in top-right corner
    if metrics_text:
        ax.text(0.98, 0.99, metrics_text, transform=ax.transAxes,
                fontsize=7, fontfamily='monospace', va='top', ha='right',
                color='#333333',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor='#CCCCCC', alpha=0.85))

    # Styling
    ax.set_xlabel("m/z", fontsize=12, fontweight='bold')
    ylabel = "Relative Intensity" if normalize else "Intensity"
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if normalize:
        # Custom y-axis labels (only show ticks at actual intensity values)
        ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
        ax.set_yticklabels(['1.0', '0.5', '0', '0.5', '1.0'])
    else:
        # Auto ticks with absolute values
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):.0f}")
        )
    ax.tick_params(axis='y', which='both', length=4)

    fig.tight_layout()

    return fig


def create_single_spectrum_plot(
    spectrum: pd.DataFrame,
    title: str = "Annotated Spectrum",
    label: str = "",
    color_map: Dict[str, str] = None,
    figsize: Tuple[int, int] = (12, 8),
    dpi: int = 300,
    show_annotations: bool = True,
    annotation_fontsize: int = 7,
    normalize: bool = True,
):
    """
    Create a single annotated spectrum plot (no mirror/reference).

    Args:
        spectrum: DataFrame with columns [mz, intensity, annotation]
        title: Plot title
        label: Label shown in top-left corner
        color_map: Custom ion color mapping
        figsize: Figure size (width, height)
        dpi: Resolution for export
        show_annotations: Whether to show ion labels
        annotation_fontsize: Font size for annotations

    Returns:
        matplotlib.figure.Figure object
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure

    if color_map is None:
        color_map = ION_COLORS

    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    # Plot peaks (positive intensities only)
    if not spectrum.empty and len(spectrum) > 0:
        df = spectrum.copy()

        if 'annotation' not in df.columns:
            df['annotation'] = 'unannotated'

        # Normalize if requested
        if normalize:
            max_int = df["intensity"].max()
            if max_int > 0:
                df["intensity"] = df["intensity"] / max_int
            ann_offset = 0.03
        else:
            ann_offset = df["intensity"].max() * 0.03 if len(df) > 0 else 1

        for _, peak in df.iterrows():
            annotation = peak.get('annotation', 'unannotated')
            color = get_ion_color(annotation, color_map)
            ax.vlines(peak['mz'], 0, peak['intensity'],
                     colors=color, linewidth=1.5, alpha=0.8)

            if show_annotations and annotation and annotation != 'unannotated':
                ax.text(peak['mz'], peak['intensity'] + ann_offset, annotation,
                       fontsize=annotation_fontsize, color=color,
                       ha='center', va='bottom', rotation=0)

    # Axis limits
    if not spectrum.empty:
        all_mz = spectrum["mz"].tolist()
        mz_min = min(all_mz) - 50
        mz_max = max(all_mz) + 50
        ax.set_xlim(mz_min, mz_max)

    if normalize:
        ax.set_ylim(0, 1.15)
    else:
        y_max = spectrum["intensity"].max() * 1.15 if not spectrum.empty else 1
        ax.set_ylim(0, y_max)

    # Label in top-left
    if label:
        ax.text(0.02, 0.99, label, transform=ax.transAxes,
                fontsize=9, fontweight='bold', va='top', ha='left',
                color='#333333',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.8))

    # Styling
    ax.set_xlabel("m/z", fontsize=12, fontweight='bold')
    ylabel = "Relative Intensity" if normalize else "Intensity"
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if normalize:
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.tick_params(axis='y', which='both', length=4)

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
