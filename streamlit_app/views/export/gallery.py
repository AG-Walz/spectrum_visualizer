"""Plot-gallery browser — called by ``export_page.py`` and ``export_actions.py`` to render the per-run expandable view of exported PNG/SVG/HTML plots, grouped by peptide variant then reference."""
from __future__ import annotations

from pathlib import Path

import streamlit as st


def _pick_browse_format(root_dir: Path) -> tuple[Path, str] | None:
    """Return ``(subdir, format)`` for inline preview — prefer SVG, fall back to PNG, ``None`` if neither exists."""
    for subdir, output_format in [(root_dir / "svg", "svg"), (root_dir / "png", "png")]:
        if subdir.exists():
            return subdir, output_format
    return None


def _group_files_by_reference(plot_files: list[Path]) -> dict[str, list[Path]]:
    """Bucket ``..._vs_<reference>.<ext>`` plot files by their reference stem."""
    ref_groups: dict[str, list[Path]] = {}
    for plot_file in plot_files:
        ref_stem = (
            plot_file.stem.rsplit("_vs_", 1)[1]
            if "_vs_" in plot_file.stem
            else "(no reference)"
        )
        ref_groups.setdefault(ref_stem, []).append(plot_file)
    return ref_groups


def _show_plot(plot_file: Path, browse_format: str) -> None:
    """Render one plot inline — ``st.image`` for PNG, raw markup for SVG."""
    if browse_format == "png":
        st.image(str(plot_file), caption=plot_file.stem, width="stretch")
    elif browse_format == "svg":
        st.markdown(plot_file.read_text(encoding="utf-8"), unsafe_allow_html=True)


def _render_download_buttons(plot_file: Path, peptide_dir_name: str, root_dir: Path) -> None:
    """Render PNG/SVG/HTML download buttons for one plot stem, skipping formats that weren't exported."""
    run_id = root_dir.name
    dl_formats = [
        ("png", "PNG", "image/png", root_dir / "png"),
        ("svg", "SVG", "image/svg+xml", root_dir / "svg"),
        ("html", "HTML", "text/html", root_dir / "html"),
    ]
    for col, (ext, label, mime, base) in zip(st.columns(3), dl_formats):
        path = base / peptide_dir_name / f"{plot_file.stem}.{ext}"
        if path.exists():
            with col:
                st.download_button(
                    label, data=path.read_bytes(),
                    file_name=f"{plot_file.stem}.{ext}", mime=mime,
                    key=f"dl_{ext}_{run_id}_{peptide_dir_name}_{plot_file.stem}",
                )


def _render_reference_group(
    ref_stem: str, ref_files: list[Path], peptide_dir_name: str,
    browse_format: str, root_dir: Path,
) -> None:
    """Render the inner expander for one reference stem, capped at 20 plots."""
    sub_label = (
        f"Reference: {ref_stem}" if ref_stem != "(no reference)"
        else "Single spectra"
    )
    with st.expander(f"{sub_label} ({len(ref_files)})", expanded=False):
        for plot_file in ref_files[:20]:
            _show_plot(plot_file, browse_format)
            _render_download_buttons(plot_file, peptide_dir_name, root_dir)
        if len(ref_files) > 20:
            st.caption(f"... and {len(ref_files) - 20} more")


def _render_peptide_dir(peptide_dir: Path, browse_format: str, root_dir: Path) -> None:
    """Render the outer expander for one peptide-variant subdirectory."""
    plot_files = sorted(peptide_dir.glob(f"*.{browse_format}"))
    if not plot_files:
        return
    with st.expander(f"{peptide_dir.name} ({len(plot_files)} plots)", expanded=False):
        for ref_stem, ref_files in sorted(_group_files_by_reference(plot_files).items()):
            _render_reference_group(ref_stem, ref_files, peptide_dir.name, browse_format, root_dir)


def _render_plot_gallery(root_dir: Path) -> None:
    """Render the gallery for plots saved under ``root_dir``, grouped peptide-variant then reference."""
    # Pick preview format: SVG if present, else PNG, else bail
    browse = _pick_browse_format(root_dir)
    if browse is None:
        return
    browse_dir, browse_format = browse

    # Discover peptide variants: one subdirectory per peptidoform
    peptide_dirs = sorted([d for d in browse_dir.iterdir() if d.is_dir()])
    if not peptide_dirs:
        return

    # Header: variant count + total plot count
    st.markdown("---")
    st.subheader("Browse Plots")
    total_plots = sum(len(list(d.glob(f"*.{browse_format}"))) for d in peptide_dirs)
    st.caption(f"{len(peptide_dirs)} peptide variants, {total_plots} plots")

    # Render one outer expander per peptide variant
    for peptide_dir in peptide_dirs:
        _render_peptide_dir(peptide_dir, browse_format, root_dir)