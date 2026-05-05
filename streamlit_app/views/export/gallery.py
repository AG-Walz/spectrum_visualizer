"""
Plot-gallery browser — per-run expandable view of exported plots.

Reads the ``png/``, ``svg/``, ``html/`` subfolders written by ``pipeline.py``
and presents them grouped by peptide / reference, each with its own download
buttons.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st


def _pick_browse_format(root_dir: Path) -> tuple[Path, str] | None:
    """Return the (dir, format) used to browse plots — prefer SVG, fall back to PNG."""
    for subdir, fmt in [(root_dir / "svg", "svg"), (root_dir / "png", "png")]:
        if subdir.exists():
            return subdir, fmt
    return None


def _group_files_by_reference(plot_files: list[Path]) -> dict[str, list[Path]]:
    """Group ``..._vs_<ref>.<ext>`` plot files by their reference stem."""
    ref_groups: dict[str, list[Path]] = {}
    for plot_file in plot_files:
        ref_stem = (
            plot_file.stem.rsplit("_vs_", 1)[1]
            if "_vs_" in plot_file.stem
            else "(no reference)"
        )
        ref_groups.setdefault(ref_stem, []).append(plot_file)
    return ref_groups


def _show_plot(plot_file: Path, browse_fmt: str) -> None:
    """Render one plot file inline (image for PNG, raw markup for SVG)."""
    if browse_fmt == "png":
        st.image(str(plot_file), caption=plot_file.stem, width="stretch")
    elif browse_fmt == "svg":
        st.markdown(plot_file.read_text(encoding="utf-8"), unsafe_allow_html=True)


def _render_download_buttons(plot_file: Path, peptide_dir_name: str, root_dir: Path) -> None:
    """Render PNG / SVG / HTML download buttons for a single plot stem."""
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
    browse_fmt: str, root_dir: Path,
) -> None:
    """Render an inner expander grouping all plots with the same reference."""
    sub_label = (
        f"Reference: {ref_stem}" if ref_stem != "(no reference)"
        else "Single spectra"
    )
    with st.expander(f"{sub_label} ({len(ref_files)})", expanded=False):
        for plot_file in ref_files[:20]:
            _show_plot(plot_file, browse_fmt)
            _render_download_buttons(plot_file, peptide_dir_name, root_dir)
        if len(ref_files) > 20:
            st.caption(f"... and {len(ref_files) - 20} more")


def _render_peptide_dir(peptide_dir: Path, browse_fmt: str, root_dir: Path) -> None:
    """Render the outer expander for one peptide variant directory."""
    plot_files = sorted(peptide_dir.glob(f"*.{browse_fmt}"))
    if not plot_files:
        return
    with st.expander(f"{peptide_dir.name} ({len(plot_files)} plots)", expanded=False):
        for ref_stem, ref_files in sorted(_group_files_by_reference(plot_files).items()):
            _render_reference_group(ref_stem, ref_files, peptide_dir.name, browse_fmt, root_dir)


def _render_plot_gallery(root_dir: Path) -> None:
    """Display plots saved in root_dir organized by peptidoform then reference."""
    browse = _pick_browse_format(root_dir)
    if browse is None:
        return
    browse_dir, browse_fmt = browse

    peptide_dirs = sorted([d for d in browse_dir.iterdir() if d.is_dir()])
    if not peptide_dirs:
        return

    st.markdown("---")
    st.subheader("Browse Plots")
    total_plots = sum(len(list(d.glob(f"*.{browse_fmt}"))) for d in peptide_dirs)
    st.caption(f"{len(peptide_dirs)} peptide variants, {total_plots} plots")

    for peptide_dir in peptide_dirs:
        _render_peptide_dir(peptide_dir, browse_fmt, root_dir)