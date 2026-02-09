"""
Export functionality for mirror plots.

Exports plots to PNG and SVG with organized folder structure.
"""
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from io import BytesIO

logger = logging.getLogger(__name__)


def sanitize_filename(name: str) -> str:
    """
    Sanitize a string for use as a filename.

    Removes or replaces invalid characters.
    """
    # Replace problematic characters
    sanitized = re.sub(r'[<>:"/\\|?*\[\]]', '_', name)
    # Replace multiple underscores with single
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores and spaces
    sanitized = sanitized.strip('_ ')
    # Truncate if too long
    if len(sanitized) > 100:
        sanitized = sanitized[:100]
    return sanitized


def create_export_directory(
    base_dir: str,
    experiment_name: str = None
) -> Tuple[Path, Path, Path]:
    """
    Create export directory structure.

    Args:
        base_dir: Base output directory
        experiment_name: Optional name for the experiment

    Returns:
        Tuple of (root_dir, png_dir, svg_dir)
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if experiment_name:
        folder_name = f"MSpecViz_{sanitize_filename(experiment_name)}_{timestamp}"
    else:
        folder_name = f"MSpecViz_export_{timestamp}"

    root_dir = Path(base_dir) / folder_name
    png_dir = root_dir / "png"
    svg_dir = root_dir / "svg"

    root_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(exist_ok=True)
    svg_dir.mkdir(exist_ok=True)

    return root_dir, png_dir, svg_dir


def export_figure(
    fig,
    output_path: Path,
    format: str = 'png',
    dpi: int = 300
) -> bool:
    """
    Export a matplotlib figure to file.

    Args:
        fig: matplotlib Figure object
        output_path: Output file path
        format: 'png' or 'svg'
        dpi: Resolution for PNG

    Returns:
        True if successful
    """
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)

        kwargs = {'format': format, 'bbox_inches': 'tight', 'facecolor': 'white'}
        if format == 'png':
            kwargs['dpi'] = dpi

        fig.savefig(str(output_path), **kwargs)
        logger.debug(f"Exported {format.upper()}: {output_path}")
        return True

    except Exception as e:
        logger.error(f"Failed to export figure: {e}")
        return False


def export_mirror_plots(
    variant_groups: Dict,
    output_dir: str,
    formats: List[str] = None,
    ion_types: List[str] = None,
    mz_tolerance: float = 0.02,
    color_map: Dict[str, str] = None,
    dpi: int = 300,
    figsize: Tuple[int, int] = (12, 8),
    progress_callback: Optional[callable] = None
) -> Dict[str, List[str]]:
    """
    Export mirror plots for all variant groups.

    Creates organized folder structure:
    output_dir/
    ├── png/
    │   ├── PEPTIDE/
    │   │   └── PEPTIDE_query1_vs_ref1.png
    │   └── PEPTIDE_Oxidation/
    │       └── ...
    └── svg/
        └── (same structure)

    Args:
        variant_groups: Dict of variant_key -> VariantGroup with loaded spectra
        output_dir: Base output directory
        formats: List of formats to export ('png', 'svg')
        ion_types: Ion types for annotation
        mz_tolerance: Tolerance for peak matching
        color_map: Custom ion colors
        dpi: PNG resolution
        figsize: Figure size
        progress_callback: Optional callback(current, total, message)

    Returns:
        Dict mapping format -> list of exported file paths
    """
    from src.peak_annotation import prepare_spectrum_for_plot
    from src.mirror_plot import create_mirror_plot, ION_COLORS

    if formats is None:
        formats = ['png', 'svg']

    if color_map is None:
        color_map = ION_COLORS

    if ion_types is None:
        ion_types = ['b', 'y']

    # Create directory structure
    root_dir, png_dir, svg_dir = create_export_directory(output_dir)

    exported_files = {fmt: [] for fmt in formats}

    # Count total plots to generate
    total_plots = 0
    for vg in variant_groups.values():
        n_query = len(vg.query_spectra)
        n_ref = len(vg.reference_spectra)
        if n_query > 0 and n_ref > 0:
            total_plots += n_query * n_ref
        elif n_query >= 2:
            # Compare query spectra against each other
            total_plots += n_query - 1

    current = 0

    # Process each variant group
    for variant_key, vg in variant_groups.items():
        # Create variant folder name
        if vg.modification_key:
            folder_name = f"{vg.base_sequence}_{sanitize_filename(vg.modification_key)}"
        else:
            folder_name = vg.base_sequence

        query_spectra = vg.query_spectra
        ref_spectra = vg.reference_spectra

        # Determine comparison pairs
        pairs = []
        if query_spectra and ref_spectra:
            # Query vs Reference
            for q_spec in query_spectra:
                for r_spec in ref_spectra:
                    pairs.append((q_spec, r_spec, "query_vs_ref"))
        elif len(query_spectra) >= 2:
            # Compare first query against rest
            first = query_spectra[0]
            for other in query_spectra[1:]:
                pairs.append((other, first, "query_vs_query"))

        # Generate plots for each pair
        for upper_spec, lower_spec, comparison_type in pairs:
            if progress_callback:
                progress_callback(current, total_plots,
                    f"Generating {folder_name}: {upper_spec.source_file} vs {lower_spec.source_file}")

            try:
                # Prepare spectra for plotting
                upper_df = prepare_spectrum_for_plot(
                    upper_spec.mz, upper_spec.intensity,
                    upper_spec.peptidoform, upper_spec.precursor_charge,
                    ion_types=ion_types, mz_tolerance=mz_tolerance
                )

                lower_df = prepare_spectrum_for_plot(
                    lower_spec.mz, lower_spec.intensity,
                    lower_spec.peptidoform, lower_spec.precursor_charge,
                    ion_types=ion_types, mz_tolerance=mz_tolerance
                )

                # Create title
                title = f"{vg.display_label}\n{Path(upper_spec.source_file).stem} vs {Path(lower_spec.source_file).stem}"

                # Create labels
                if comparison_type == "query_vs_ref":
                    upper_label = f"Query: {Path(upper_spec.source_file).stem}"
                    lower_label = f"Reference: {Path(lower_spec.source_file).stem}"
                else:
                    upper_label = f"{Path(upper_spec.source_file).stem}"
                    lower_label = f"{Path(lower_spec.source_file).stem}"

                # Create figure
                fig = create_mirror_plot(
                    upper_df, lower_df,
                    title=title,
                    upper_label=upper_label,
                    lower_label=lower_label,
                    color_map=color_map,
                    figsize=figsize,
                    dpi=dpi
                )

                # Generate filename
                upper_name = sanitize_filename(Path(upper_spec.source_file).stem)
                lower_name = sanitize_filename(Path(lower_spec.source_file).stem)
                # Truncate if needed
                if len(upper_name) > 30:
                    upper_name = upper_name[:30]
                if len(lower_name) > 30:
                    lower_name = lower_name[:30]

                base_filename = f"{vg.base_sequence}_{upper_name}_vs_{lower_name}"

                # Export to each format
                for fmt in formats:
                    if fmt == 'png':
                        dir_path = png_dir / folder_name
                    else:
                        dir_path = svg_dir / folder_name

                    file_path = dir_path / f"{base_filename}.{fmt}"

                    if export_figure(fig, file_path, format=fmt, dpi=dpi):
                        exported_files[fmt].append(str(file_path))

                # Close figure to free memory
                import matplotlib.pyplot as plt
                plt.close(fig)

            except Exception as e:
                logger.error(f"Error generating plot for {variant_key}: {e}")

            current += 1

    if progress_callback:
        progress_callback(total_plots, total_plots, "Export complete")

    # Log summary
    for fmt, files in exported_files.items():
        logger.info(f"Exported {len(files)} {fmt.upper()} files to {root_dir}")

    return exported_files


def create_summary_excel(
    variant_groups: Dict,
    output_path: str,
    include_stats: bool = True
) -> bool:
    """
    Create an Excel summary of the analysis.

    Args:
        variant_groups: Dict of variant_key -> VariantGroup
        output_path: Path for Excel file
        include_stats: Whether to include annotation statistics

    Returns:
        True if successful
    """
    try:
        import pandas as pd

        rows = []
        for variant_key, vg in variant_groups.items():
            base_row = {
                'Peptide': vg.base_sequence,
                'Modification': vg.modification_key or "(none)",
                'Display Label': vg.display_label,
                'Total PSMs': len(vg.matches),
                'Query Spectra': len(vg.query_spectra),
                'Reference Spectra': len(vg.reference_spectra)
            }

            # Add per-spectrum info
            for i, spec in enumerate(vg.query_spectra):
                row = base_row.copy()
                row['Spectrum Type'] = 'Query'
                row['Spectrum Index'] = i + 1
                row['Source File'] = spec.source_file
                row['Spectrum ID'] = spec.spectrum_id
                row['Precursor m/z'] = spec.precursor_mz
                row['Charge'] = spec.precursor_charge
                row['RT'] = spec.rt
                rows.append(row)

            for i, spec in enumerate(vg.reference_spectra):
                row = base_row.copy()
                row['Spectrum Type'] = 'Reference'
                row['Spectrum Index'] = i + 1
                row['Source File'] = spec.source_file
                row['Spectrum ID'] = spec.spectrum_id
                row['Precursor m/z'] = spec.precursor_mz
                row['Charge'] = spec.precursor_charge
                row['RT'] = spec.rt
                rows.append(row)

        df = pd.DataFrame(rows)
        df.to_excel(output_path, index=False)
        logger.info(f"Saved summary Excel: {output_path}")
        return True

    except Exception as e:
        logger.error(f"Failed to create Excel summary: {e}")
        return False
