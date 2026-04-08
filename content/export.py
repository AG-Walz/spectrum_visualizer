"""
Export Page — Parquet, CSV/XLSX, mirror plots, summary plots, HTML, ZIP download.

Export controls, history, gallery, and ZIP download are each @st.fragment
so interactions in one section don't rerun the others.
"""
import io
import json
import logging
import shutil
import zipfile
import streamlit as st
from pathlib import Path
from datetime import datetime

from src.common.common import page_setup, save_params
from content._shared import (
    STATE_PREFIX, get_state_key, restore_analysis_if_needed,
    build_run_settings, apply_run_settings,
)
from utils.helpers import parse_peptide_list

params = page_setup()
logger = logging.getLogger(__name__)


def _p(key, default=None):
    """Get a param value: session_state first (if widget present), then params.json."""
    return st.session_state.get(key, params.get(key, default))


# ---------------------------------------------------------------------------
# Mirror plot export helper
# ---------------------------------------------------------------------------

def _export_mirror_plots(export_root: Path, do_html: bool):
    """Generate mirror/single plots organized by peptidoform.

    Returns (file_count, list_of_file_paths).
    """
    from src.spectrum_metrics import compute_all_metrics, format_metrics_text
    from src.summary_plots import sanitize_filename
    from src.mirror_plot import (
        create_mirror_plot, create_single_spectrum_plot,
        figure_to_bytes, create_plotly_mirror_plot,
    )
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    annotated_queries = st.session_state.get(get_state_key('annotated_queries'), {})
    reference_store = st.session_state.get(get_state_key('reference_store'), {})
    variant_groups = st.session_state.get(get_state_key('variant_groups'), {})

    if not annotated_queries:
        st.warning("No annotated spectra in session -- re-run **Search & Load** first.")
        return 0, []

    st.caption(f"  {len(annotated_queries)} annotated spectra in session")

    normalize = _p('intensity-mode', 'Relative') == 'Relative'

    png_dir = export_root / "png"
    svg_dir = export_root / "svg"
    html_dir = export_root / "html" if do_html else None

    from src.spectrum_loader import spectrum_lookup_key
    psm_id_to_vk = {}
    for vk, vg in variant_groups.items():
        for match in vg.matches:
            key = f"{match.spectrum_id}|{match.peptidoform}"
            for pid, (spec, _) in annotated_queries.items():
                if spectrum_lookup_key(spec) == key:
                    psm_id_to_vk[pid] = vk
                    break

    saved_files = []
    errors = []

    for pid, (query_spec, query_df) in sorted(annotated_queries.items()):
        vk = psm_id_to_vk.get(pid)
        refs = reference_store.get(vk, []) if vk else []
        folder_name = sanitize_filename(query_spec.peptidoform)

        if refs:
            for ref_spec, ref_df in refs:
                try:
                    metrics = compute_all_metrics(query_df, ref_df)
                    metrics_text = format_metrics_text(metrics)
                    ref_source = Path(ref_spec.source_file).stem if ref_spec.source_file else "Koina"
                    safe_name = sanitize_filename(
                        f"{Path(query_spec.source_file).stem}_vs_{ref_source}"
                    )
                    title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
                    query_label = f"Query: {Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"
                    ref_label = f"Ref: {ref_source} (RT: {ref_spec.rt:.1f})"

                    fig = create_mirror_plot(
                        upper_spectrum=query_df, lower_spectrum=ref_df,
                        title=title, upper_label=query_label, lower_label=ref_label,
                        metrics_text=metrics_text, normalize=normalize,
                    )

                    for fmt_dir, fmt in [(png_dir, "png"), (svg_dir, "svg")]:
                        out_dir = fmt_dir / folder_name
                        out_dir.mkdir(parents=True, exist_ok=True)
                        out_file = out_dir / f"{safe_name}.{fmt}"
                        out_file.write_bytes(figure_to_bytes(fig, format=fmt))
                        saved_files.append(str(out_file))

                    plt.close(fig)

                    if do_html and html_dir:
                        plotly_fig = create_plotly_mirror_plot(
                            upper_spectrum=query_df, lower_spectrum=ref_df,
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
                    logger.error(f"Mirror plot failed for PSM {pid}: {e}")
        else:
            try:
                safe_name = sanitize_filename(Path(query_spec.source_file).stem)
                title = f"{query_spec.peptidoform} +{query_spec.precursor_charge}"
                label = f"{Path(query_spec.source_file).stem} (RT: {query_spec.rt:.1f})"

                fig = create_single_spectrum_plot(
                    spectrum=query_df, title=title, label=label, normalize=normalize,
                )

                for fmt_dir, fmt in [(png_dir, "png"), (svg_dir, "svg")]:
                    out_dir = fmt_dir / folder_name
                    out_dir.mkdir(parents=True, exist_ok=True)
                    out_file = out_dir / f"{safe_name}.{fmt}"
                    out_file.write_bytes(figure_to_bytes(fig, format=fmt))
                    saved_files.append(str(out_file))

                plt.close(fig)

            except Exception as e:
                errors.append(f"{folder_name}/{safe_name}: {e}")
                logger.error(f"Single plot failed for PSM {pid}: {e}")

    if errors:
        st.warning(f"{len(errors)} plot(s) failed:")
        for err in errors[:10]:
            st.caption(f"  {err}")

    return len(saved_files), saved_files


# ---------------------------------------------------------------------------
# Export execution
# ---------------------------------------------------------------------------

def _run_export(export_dir: str, do_parquet: bool, do_csv: bool,
                do_mirror: bool, do_summary: bool, do_html: bool):
    """Execute all selected exports with versioned directory."""
    import pandas as pd
    import traceback as _tb

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    export_root = Path(export_dir) / f"export_{timestamp}"
    export_root.mkdir(parents=True, exist_ok=True)

    saved_files = []

    with st.status("Exporting...", expanded=True) as status:
        run_settings = build_run_settings()
        settings_path = export_root / "run_settings.json"
        with open(settings_path, "w") as f:
            json.dump(run_settings, f, indent=2, default=str)
        saved_files.append(str(settings_path))

        if do_parquet:
            st.caption("Saving Parquet DataFrames...")
            for name in ('psm_lf', 'peaks_lf', 'seq_lf', 'heatmap_lf'):
                lf = st.session_state.get(get_state_key(name))
                if lf is not None:
                    out_file = export_root / f"{name.replace('_lf', '')}.parquet"
                    lf.collect().write_parquet(out_file)
                    saved_files.append(str(out_file))

        if do_csv:
            comparison_df = st.session_state.get(get_state_key('comparison_df'))
            if comparison_df is not None and not comparison_df.empty:
                st.caption("Saving comparison summary...")
                csv_path = export_root / "comparison_summary.csv"
                comparison_df.to_csv(csv_path, index=False)
                saved_files.append(str(csv_path))

                try:
                    xlsx_path = export_root / "comparison_summary.xlsx"
                    comparison_df.to_excel(xlsx_path, index=False, engine="openpyxl")
                    saved_files.append(str(xlsx_path))
                except ImportError:
                    logger.warning("openpyxl not installed -- skipping XLSX export")
                except Exception as e:
                    logger.warning(f"XLSX export failed: {e}")

        if do_mirror or do_html:
            st.caption("Generating mirror plots...")
            try:
                n_plots, plot_files = _export_mirror_plots(export_root, do_html)
                saved_files.extend(plot_files)
                st.caption(f"  Generated {n_plots} mirror plot file(s)")
            except Exception as e:
                st.error(f"Mirror plot export failed: {e}")
                logger.error(f"Mirror plot export failed:\n{_tb.format_exc()}")
        else:
            st.caption("Mirror plots: skipped (unchecked)")

        if do_summary:
            comparison_df = st.session_state.get(get_state_key('comparison_df'))
            if comparison_df is not None and not comparison_df.empty:
                st.caption("Generating summary plots...")
                from src.summary_plots import generate_all_summary_plots
                try:
                    results = generate_all_summary_plots(
                        source=comparison_df,
                        output_dir=str(export_root),
                        formats=["png", "svg"],
                    )
                    n_summary = sum(len(v) for v in results.values())
                    st.caption(f"  Generated {n_summary} summary plot(s)")
                    for pt, paths in results.items():
                        for fmt, p in paths.items():
                            saved_files.append(p)
                except Exception as e:
                    logger.error(f"Summary plot generation failed: {e}")
                    st.warning(f"Summary plots failed: {e}")

        status.update(label=f"Exported {len(saved_files)} files to {export_root}", state="complete")

    if saved_files:
        with st.expander(f"Exported files ({len(saved_files)})", expanded=False):
            for f in saved_files:
                st.caption(f"  {f}")

    _render_plot_gallery(export_root)


# ---------------------------------------------------------------------------
# Plot gallery
# ---------------------------------------------------------------------------

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
            for pf in plot_files:
                ref_stem = pf.stem.rsplit("_vs_", 1)[1] if "_vs_" in pf.stem else "(no reference)"
                ref_groups.setdefault(ref_stem, []).append(pf)

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
                        dl_col1, dl_col2, dl_col3 = st.columns(3)
                        with dl_col1:
                            png_path = root_dir / "png" / peptide_dir.name / f"{plot_file.stem}.png"
                            if png_path.exists():
                                st.download_button(
                                    "PNG", data=png_path.read_bytes(),
                                    file_name=f"{plot_file.stem}.png",
                                    mime="image/png",
                                    key=f"dl_png_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )
                        with dl_col2:
                            svg_path = root_dir / "svg" / peptide_dir.name / f"{plot_file.stem}.svg"
                            if svg_path.exists():
                                st.download_button(
                                    "SVG", data=svg_path.read_bytes(),
                                    file_name=f"{plot_file.stem}.svg",
                                    mime="image/svg+xml",
                                    key=f"dl_svg_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )
                        with dl_col3:
                            html_file = html_base_dir / peptide_dir.name / f"{plot_file.stem}.html"
                            if html_file.exists():
                                st.download_button(
                                    "HTML", data=html_file.read_bytes(),
                                    file_name=f"{plot_file.stem}.html",
                                    mime="text/html",
                                    key=f"dl_html_{run_id}_{peptide_dir.name}_{plot_file.stem}",
                                )

                    if len(ref_files) > 20:
                        st.caption(f"... and {len(ref_files) - 20} more")


# ---------------------------------------------------------------------------
# Export history
# ---------------------------------------------------------------------------

@st.fragment
def export_history_fragment(export_dir: str):
    """Browse previous export runs, reload their settings, or delete them."""
    st.subheader("Export History")

    export_path = Path(export_dir)
    if not export_path.exists():
        st.caption(f"No exports directory found at `{export_dir}`.")
        return

    runs = []
    for settings_file in sorted(export_path.glob("*/run_settings.json"), reverse=True):
        try:
            with open(settings_file) as f:
                run_cfg = json.load(f)
            runs.append((settings_file.parent, run_cfg))
        except Exception:
            continue

    if not runs:
        st.caption(f"No completed export runs found in `{export_dir}`.")
        return

    st.caption(f"{len(runs)} run(s) in `{export_dir}`")

    for run_dir, run_cfg in runs:
        peptides = run_cfg.get('peptides', [])
        idxml_count = len(run_cfg.get('idxml_files', []))
        mzml_count = len(run_cfg.get('mzml_files', []))
        created = (run_cfg.get('_created_at') or '')[:10]

        pep_preview = ', '.join(peptides[:3])
        if len(peptides) > 3:
            pep_preview += f' ... +{len(peptides) - 3}'

        label = (
            f"**{run_dir.name}**"
            f" -- {created}"
            f" -- {len(peptides)} peptide(s)"
            f" -- {idxml_count} idXML / {mzml_count} mzML"
        )

        with st.expander(label, expanded=False):
            col_meta, col_actions = st.columns([3, 1])

            with col_meta:
                st.markdown(f"**Peptides:** {pep_preview or '--'}")
                st.markdown(f"**Files:** {idxml_count} idXML, {mzml_count} mzML")
                ref_sub = run_cfg.get('reference_substring') or '--'
                st.markdown(f"**Reference substring:** `{ref_sub}`")
                st.markdown(f"**Ion types:** {', '.join(run_cfg.get('ion_types', []))}")
                st.markdown(f"**m/z tolerance:** {run_cfg.get('mz_tolerance', '--')} Da")
                created_full = (run_cfg.get('_created_at') or '--')[:19]
                last_run_full = (run_cfg.get('_last_run') or '--')[:19]
                st.caption(f"Created: {created_full}  |  Last run: {last_run_full}")

            with col_actions:
                if st.button("Reload settings", key=f"reload_{run_dir.name}"):
                    apply_run_settings(run_cfg)
                    st.success("Settings loaded -- switch to Setup page to review.")
                    st.rerun()

                confirm_key = f"confirm_delete_{run_dir.name}"
                if st.session_state.get(confirm_key):
                    st.warning("Delete all files in this run?")
                    c1, c2 = st.columns(2)
                    with c1:
                        if st.button("Yes, delete", key=f"yes_del_{run_dir.name}", type="primary"):
                            shutil.rmtree(run_dir)
                            del st.session_state[confirm_key]
                            st.rerun()
                    with c2:
                        if st.button("Cancel", key=f"cancel_del_{run_dir.name}"):
                            del st.session_state[confirm_key]
                            st.rerun()
                else:
                    if st.button("Delete run", key=f"delete_{run_dir.name}"):
                        st.session_state[confirm_key] = True
                        st.rerun()

            _render_plot_gallery(run_dir)


# ---------------------------------------------------------------------------
# ZIP download
# ---------------------------------------------------------------------------

@st.fragment
def zip_download_fragment(export_root: Path) -> None:
    """Deferred zip download of the export directory."""
    def _make_zip() -> bytes:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for file in export_root.rglob("*"):
                if not file.is_file():
                    continue
                zf.write(file, file.relative_to(export_root.parent))
        return buf.getvalue()

    st.download_button(
        label="Download Export as ZIP",
        data=_make_zip,
        file_name=f"{export_root.name}.zip",
        mime="application/zip",
    )


# ---------------------------------------------------------------------------
# Export controls fragment
# ---------------------------------------------------------------------------

@st.fragment
def export_controls_fragment():
    """Export format selection and trigger button."""
    if not st.session_state.get(get_state_key('data_ready')):
        st.info("Load data first (Search & Load page).")
        return

    _ws = st.session_state.get('workspace', None)
    _export_fallback = str(Path(_ws) / "exports") if _ws else "./exports"
    export_default = params.get("export-dir", "") or _export_fallback
    st.text_input(
        "Export Directory",
        value=export_default,
        key="export-dir"
    )

    export_col1, export_col2 = st.columns(2)
    with export_col1:
        export_parquet = st.checkbox("Parquet DataFrames", value=True,
                                     key=f"{STATE_PREFIX}_exp_parquet")
        export_csv = st.checkbox("Comparison Summary (CSV/XLSX)", value=True,
                                 key=f"{STATE_PREFIX}_exp_csv")
        export_html = st.checkbox("Interactive Plots (HTML)", value=False,
                                  key=f"{STATE_PREFIX}_exp_html")
    with export_col2:
        export_mirror = st.checkbox("Mirror Plots (PNG/SVG)", value=True,
                                    key=f"{STATE_PREFIX}_exp_mirror")
        export_summary = st.checkbox("Summary Analysis Plots", value=True,
                                     key=f"{STATE_PREFIX}_exp_summary")

    if st.button("Export All", key=f"{STATE_PREFIX}_export_btn", type="primary"):
        _run_export(export_dir, export_parquet, export_csv, export_mirror,
                    export_summary, export_html)


# ---------------------------------------------------------------------------
# Page render
# ---------------------------------------------------------------------------

def render():
    restore_analysis_if_needed()
    st.title("Export")
    st.caption("Export analysis results, plots, and data files")

    export_controls_fragment()

    _ws = st.session_state.get('workspace', None)
    _fallback = str(Path(_ws) / "exports") if _ws else "./exports"
    export_dir = _p('export-dir', '') or _fallback

    st.markdown("---")
    export_history_fragment(export_dir)

    # ZIP download
    export_path = Path(export_dir)
    if export_path.exists() and any(export_path.rglob("*")):
        st.markdown("---")
        st.subheader("Download")
        zip_download_fragment(export_path)

    save_params(params)


render()
