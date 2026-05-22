"""Summary-plot orchestrator: aggregates the comparison table by peptidoform and dispatches to summary_dual_axis / summary_ridgeline.
Role — top-level coordinator for the summary-plot layer; the only module here
that touches both the raw comparison data and the per-type renderers. What it
produces — prepare_summary_data turns a comparison table into a per-peptidoform
aggregate frame (count, median, mean, std, 95% confidence interval, dominant
charge); generate_all_summary_plots runs that aggregation once and writes the
requested figure types to <output_dir>/summary/ in each requested file format.
Call relationships — generate_all_summary_plots is called by
core.export.export_pipeline (and indirectly cli.main); prepare_summary_data
is also called directly by streamlit_app.views.explore.explore_page for the
in-app preview. Dispatches to plot_classic_dual_axis and plot_ridgeline,
and draws metric names/labels from core.plotting._summary_shared.
"""
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from core.defaults import DEFAULTS
from core.plotting._summary_shared import METRIC_COLUMNS, _metric_label
from core.plotting.summary_dual_axis import plot_classic_dual_axis
from core.plotting.summary_ridgeline import plot_ridgeline

logger = logging.getLogger(__name__)

def prepare_summary_data(
    source: str | Path | pd.DataFrame,
    metric: str = DEFAULTS.summary_metric,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the comparison table and aggregate per-peptidoform stats for metric.
    source is either a CSV path or an already-built comparison DataFrame;
    the latter is copied so callers' frames are never mutated. A peptidoform
    label is synthesised per row from the modification string (falling back to
    the bare Peptide when unmodified) joined with the precursor Charge.
    Rows with a missing value in metric are dropped before aggregation.
    Returns (raw_frame, aggregated_frame): raw_frame is the per-row table
    augmented with the peptidoform column, and aggregated_frame carries
    one row per peptidoform with count, median, mean, std, q25,
    q75, the 95% confidence interval (ci_lo/ci_hi), and the dominant
    charge — sorted ascending by count.
    Raises ValueError when metric is not a column of the input table.
    """
    # Load: accept either an on-disk CSV path or an already-built DataFrame
    if isinstance(source, (str, Path)):
        raw_frame = pd.read_csv(source, low_memory=False)
    else:
        raw_frame = source.copy()

    # Validate: fail fast if the requested metric column is missing
    if metric not in raw_frame.columns:
        raise ValueError(
            f"Metric '{metric}' not found. Available: "
            f"{[c for c in raw_frame.columns if c in METRIC_COLUMNS]}"
        )

    # Build peptidoform label: modification-or-base sequence joined with charge
    modification = raw_frame["Modification"].fillna("").astype(str)
    inline = modification.mask(
        modification.isin(["", "(unmodified)"]), raw_frame["Peptide"].astype(str),
    )
    raw_frame["peptidoform"] = inline + "/" + raw_frame["Charge"].astype(str)

    raw_frame = raw_frame.dropna(subset=[metric]).copy()

    # Aggregate: per-peptidoform count/quartile stats + 95% CI + dominant charge
    grouped = raw_frame.groupby("peptidoform", sort=False)
    agg = grouped[metric].agg([
        "count", "median", "mean", "std",
        ("q25", lambda x: x.quantile(0.25)),
        ("q75", lambda x: x.quantile(0.75)),
    ])
    agg["std"] = agg["std"].fillna(0)
    agg["ci_lo"] = agg["mean"] - 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["ci_hi"] = agg["mean"] + 1.96 * agg["std"] / np.sqrt(agg["count"])
    agg["charge"] = grouped["Charge"].agg(lambda x: x.mode().iloc[0])

    # Finalise: sort ascending by count so plots draw least→most prevalent
    agg = agg.reset_index().sort_values("count", ascending=True).reset_index(drop=True)
    return raw_frame, agg


def generate_all_summary_plots(
    source: str | Path | pd.DataFrame,
    output_dir: str | Path,
    metric: str = DEFAULTS.summary_metric,
    plot_types: list[str] | None = None,
    formats: list[str] | None = None,
    export_resolution: int = 300,
) -> dict[str, dict[str, str]]:
    """Render every requested summary plot type into <output_dir>/summary/.
    source is passed straight through to prepare_summary_data, which runs
    once so the aggregation is shared across all plot types. plot_types defaults
    to DEFAULTS.summary_plot_types and formats defaults to ["png"];
    each recognised type (dual_axis, ridgeline) is rendered once and saved
    under every requested format. The ridgeline figure for the default metric is
    given the special basename ridgeline_r2_binned.
    Unknown plot types, renderer ValueErrors, and per-file OSErrors are
    logged and skipped rather than raised, so one bad type/format does not abort
    the rest. Returns {plot_type: {format: path}} covering only the files
    that were written successfully.
    """
    # Resolve defaults: fall back to configured plot types + PNG-only output
    if plot_types is None:
        plot_types = list(DEFAULTS.summary_plot_types)
    if formats is None:
        formats = ["png"]

    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Prep data: aggregate once, reuse across every requested plot type
    raw_frame, aggregated_frame = prepare_summary_data(source, metric=metric)
    logger.info(f"Summary data: {len(aggregated_frame)} peptidoforms, {len(raw_frame)} rows, metric={metric}")

    # Render & save: dispatch each plot type directly, write per-format files
    results: dict[str, dict[str, str]] = {}
    for plot_type in plot_types:
        if plot_type == "dual_axis":
            renderer = plot_classic_dual_axis
        elif plot_type == "ridgeline":
            renderer = plot_ridgeline
        else:
            logger.warning("Unknown plot type: %s (skipping)", plot_type)
            continue

        try:
            fig = renderer(raw_frame, aggregated_frame, metric=metric,
                           export_resolution=export_resolution, title=None)
        except ValueError as error:
            logger.error("Failed to create %s plot: %s", plot_type, error)
            continue

        basename = ("ridgeline_r2_binned"
                    if plot_type == "ridgeline" and metric == DEFAULTS.summary_metric
                    else plot_type)
        paths: dict[str, str] = {}
        for output_format in formats:
            output_path = summary_dir / f"{basename}.{output_format}"
            try:
                fig.savefig(str(output_path), format=output_format, dpi=export_resolution,
                            bbox_inches="tight", facecolor="white")
                paths[output_format] = str(output_path)
            except OSError as error:
                logger.error("Failed to save %s.%s: %s", basename, output_format, error)
        plt.close(fig)
        if paths:
            results[plot_type] = paths
    return results
