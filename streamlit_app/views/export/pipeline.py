"""
Export pipeline — actually writes files to disk.

``_run_export`` is the orchestration entry point called from the Export button.
It writes ``run_settings.json``, Parquet frames, CSV/XLSX, mirror plots, and
summary plots into a fresh ``export_<timestamp>/`` directory.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

from streamlit_app.utils.state import get_param, get_state_key
from streamlit_app.utils.page_state import build_run_settings
from streamlit_app.views.export.gallery import _render_plot_gallery
from mspecviz.annotation import compute_all_metrics, format_metrics_text
from mspecviz.common.defaults import DEFAULTS
from mspecviz.common.text_utils import sanitize_filename
from mspecviz.pipeline import (
    build_psm_to_variant_mapping,
    format_spectrum_labels,
    get_source_label,
)
from mspecviz.plotting import (
    create_mirror_plot,
    create_plotly_mirror_plot,
    create_single_spectrum_plot,
    figure_to_bytes,
    generate_all_summary_plots,
)

logger = logging.getLogger(__name__)


def _save_figure(fig, png_dir, svg_dir, folder_name, safe_name, saved_files):
    """Save a matplotlib figure to PNG and SVG, then close it."""
    for fmt_dir, fmt in [(png_dir, "png"), (svg_dir, "svg")]:
        out_dir = fmt_dir / folder_name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"{safe_name}.{fmt}"
        out_file.write_bytes(figure_to_bytes(fig, format=fmt))
        saved_files.append(str(out_file))
    plt.close(fig)


def _export_spectrum_plots(params: dict, export_root: Path, do_html: bool):
    """Generate mirror/single plots organized by peptidoform.

    Returns (file_count, list_of_file_paths).
    """
    annotated_queries = st.session_state.get(get_state_key("annotated_queries"), {})
    reference_store = st.session_state.get(get_state_key("reference_store"), {})
    variant_groups = st.session_state.get(get_state_key("variant_groups"), {})

    if not annotated_queries:
        st.warning("No annotated spectra in session -- re-run **Search & Load** first.")
        return 0, []

    st.caption(f"  {len(annotated_queries)} annotated spectra in session")

    normalize = get_param(params, "intensity_mode", DEFAULTS.intensity_mode) == "Relative"
    png_dir = export_root / "png"
    svg_dir = export_root / "svg"
    html_dir = export_root / "html" if do_html else None

    psm_id_to_variant_key = build_psm_to_variant_mapping(variant_groups, annotated_queries)

    saved_files: list[str] = []
    errors: list[str] = []

    for psm_id, (query_spectrum, query_frame) in sorted(annotated_queries.items()):
        variant_key = psm_id_to_variant_key.get(psm_id)
        references = reference_store.get(variant_key, []) if variant_key else []
        folder_name = sanitize_filename(query_spectrum.peptidoform)
        title = f"{query_spectrum.peptidoform} +{query_spectrum.precursor_charge}"

        if references:
            for reference_spectrum, reference_frame in references:
                try:
                    metrics = compute_all_metrics(query_frame, reference_frame)
                    safe_name = sanitize_filename(
                        f"{get_source_label(query_spectrum)}_vs_{get_source_label(reference_spectrum)}"
                    )
                    query_label, ref_label = format_spectrum_labels(query_spectrum, reference_spectrum)
                    metrics_text = format_metrics_text(metrics)

                    fig = create_mirror_plot(
                        upper_spectrum=query_frame, lower_spectrum=reference_frame,
                        title=title, upper_label=query_label, lower_label=ref_label,
                        metrics_text=metrics_text, normalize=normalize,
                    )
                    _save_figure(fig, png_dir, svg_dir, folder_name, safe_name, saved_files)

                    if do_html and html_dir:
                        plotly_fig = create_plotly_mirror_plot(
                            upper_spectrum=query_frame, lower_spectrum=reference_frame,
                            title=title, upper_label=query_label,
                            lower_label=ref_label, metrics_text=metrics_text,
                            normalize=normalize,
                        )
                        out_dir = html_dir / folder_name
                        out_dir.mkdir(parents=True, exist_ok=True)
                        html_path = out_dir / f"{safe_name}.html"
                        plotly_fig.write_html(str(html_path), include_plotlyjs="cdn")
                        saved_files.append(str(html_path))

                except Exception as e:
                    errors.append(f"{folder_name}/{safe_name}: {e}")
                    logger.error(f"Mirror plot failed for PSM {psm_id}: {e}")
        else:
            try:
                safe_name = sanitize_filename(Path(query_spectrum.source_file).stem)
                label = f"{Path(query_spectrum.source_file).stem} (RT: {query_spectrum.rt:.1f})"
                fig = create_single_spectrum_plot(
                    spectrum=query_frame, title=title, label=label, normalize=normalize,
                )
                _save_figure(fig, png_dir, svg_dir, folder_name, safe_name, saved_files)

            except Exception as e:
                errors.append(f"{folder_name}/{safe_name}: {e}")
                logger.error(f"Single plot failed for PSM {psm_id}: {e}")

    if errors:
        st.warning(f"{len(errors)} plot(s) failed:")
        for err in errors[:10]:
            st.caption(f"  {err}")

    return len(saved_files), saved_files


def _run_export(params: dict, export_dir: str, do_parquet: bool, do_csv: bool,
                do_mirror: bool, do_summary: bool, do_html: bool) -> Path:
    """Execute all selected exports with versioned directory.

    Returns the path of the newly created export root so the caller can then
    render the plot gallery for it.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    export_root = Path(export_dir) / f"export_{timestamp}"
    export_root.mkdir(parents=True, exist_ok=True)

    saved_files: list[str] = []

    with st.status("Exporting...", expanded=True) as status:
        run_settings = build_run_settings()
        settings_path = export_root / "run_settings.json"
        with open(settings_path, "w") as f:
            json.dump(run_settings, f, indent=2, default=str)
        saved_files.append(str(settings_path))

        if do_parquet:
            st.caption("Saving Parquet DataFrames...")
            for name in ("psm_frame", "peaks_frame", "sequence_frame", "heatmap_frame"):
                frame = st.session_state.get(get_state_key(name))
                if frame is not None:
                    out_file = export_root / f"{name.replace('_frame', '')}.parquet"
                    frame.collect().write_parquet(out_file)
                    saved_files.append(str(out_file))

        if do_csv:
            comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
            if comparison_frame is not None and not comparison_frame.empty:
                st.caption("Saving comparison summary...")
                csv_path = export_root / "comparison_summary.csv"
                comparison_frame.to_csv(csv_path, index=False)
                saved_files.append(str(csv_path))

                try:
                    xlsx_path = export_root / "comparison_summary.xlsx"
                    comparison_frame.to_excel(xlsx_path, index=False, engine="openpyxl")
                    saved_files.append(str(xlsx_path))
                except ImportError:
                    logger.warning("openpyxl not installed -- skipping XLSX export")
                except Exception as e:
                    logger.warning(f"XLSX export failed: {e}")

        if do_mirror or do_html:
            st.caption("Generating mirror plots...")
            n_plots, plot_files = _export_spectrum_plots(params, export_root, do_html)
            saved_files.extend(plot_files)
            st.caption(f"  Generated {n_plots} mirror plot file(s)")
        else:
            st.caption("Mirror plots: skipped (unchecked)")

        if do_summary:
            comparison_frame = st.session_state.get(get_state_key("comparison_frame"))
            if comparison_frame is not None and not comparison_frame.empty:
                st.caption("Generating summary plots...")
                results = generate_all_summary_plots(
                    source=comparison_frame,
                    output_dir=str(export_root),
                    formats=["png", "svg"],
                )
                n_summary = sum(len(v) for v in results.values())
                st.caption(f"  Generated {n_summary} summary plot(s)")
                for pt, paths in results.items():
                    for fmt, p in paths.items():
                        saved_files.append(p)

        status.update(label=f"Exported {len(saved_files)} files to {export_root}", state="complete")

    if saved_files:
        with st.expander(f"Exported files ({len(saved_files)})", expanded=False):
            for f in saved_files:
                st.caption(f"  {f}")

    _render_plot_gallery(export_root)
    return export_root
