"""
Plot-gallery browser — per-run expandable view of exported plots.

Reads the ``png/``, ``svg/``, ``html/`` subfolders written by ``pipeline.py``
and presents them grouped by peptide / reference, each with its own download
buttons.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st


def _render_plot_gallery(root_dir: Path) -> None:
    """Display plots saved in root_dir organized by peptidoform then reference."""
    format_dirs = [(root_dir / "svg", "svg"), (root_dir / "png", "png")]
    existing_format_dirs = [(d, fmt) for d, fmt in format_dirs if d.exists()]

    if not existing_format_dirs:
        return

    browse_dir, browse_fmt = existing_format_dirs[0]
    html_base_dir = root_dir / "html"

    peptide_dirs = sorted([d for d in browse_dir.iterdir() if d.is_dir()])
    if not peptide_dirs:
        return

    st.markdown("---")
    st.subheader("Browse Plots")

    total_plots = sum(len(list(d.glob(f"*.{browse_fmt}"))) for d in peptide_dirs)
    st.caption(f"{len(peptide_dirs)} peptide variants, {total_plots} plots")

    for peptide_dir in peptide_dirs:
        plot_files = sorted(peptide_dir.glob(f"*.{browse_fmt}"))
        if not plot_files:
            continue

        with st.expander(f"{peptide_dir.name} ({len(plot_files)} plots)", expanded=False):
            ref_groups: dict = {}
            for plot_file in plot_files:
                ref_stem = (
                    plot_file.stem.rsplit("_vs_", 1)[1]
                    if "_vs_" in plot_file.stem
                    else "(no reference)"
                )
                ref_groups.setdefault(ref_stem, []).append(plot_file)

            for ref_stem, ref_files in sorted(ref_groups.items()):
                sub_label = (
                    f"Reference: {ref_stem}" if ref_stem != "(no reference)"
                    else "Single spectra"
                )
                with st.expander(f"{sub_label} ({len(ref_files)})", expanded=False):
                    for plot_file in ref_files[:20]:
                        if browse_fmt == "png":
                            st.image(str(plot_file), caption=plot_file.stem, width="stretch")
                        elif browse_fmt == "svg":
                            try:
                                svg_content = plot_file.read_text()
                                st.markdown(svg_content, unsafe_allow_html=True)
                            except Exception:
                                st.caption(f"Could not render: {plot_file.name}")

                        run_id = root_dir.name
                        dl_formats = [
                            ("png", "PNG", "image/png", root_dir / "png"),
                            ("svg", "SVG", "image/svg+xml", root_dir / "svg"),
                            ("html", "HTML", "text/html", html_base_dir),
                        ]
                        for col, (ext, label, mime, base) in zip(st.columns(3), dl_formats):
                            path = base / peptide_dir.name / f"{plot_file.stem}.{ext}"
                            if path.exists():
                                with col:
                                    st.download_button(
                                        label, data=path.read_bytes(),
                                        file_name=f"{plot_file.stem}.{ext}", mime=mime,
                                        key=f"dl_{ext}_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                    )

                    if len(ref_files) > 20:
                        st.caption(f"... and {len(ref_files) - 20} more")
