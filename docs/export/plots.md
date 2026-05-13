# Plot Outputs

What plot files an export run produces, how they are named, and the
Streamlit-only gallery/history views that let you browse past runs.

- For per-spectrum **mirror plots** see [Mirror plots](#mirror-plots).
- For per-run **summary plots** see [Summary plots](#summary-plots).
- For Streamlit's post-export browsers see [Gallery & history](#gallery--history).

For the directory tree see [layout.md](layout.md); for the
comparison-summary CSV/XLSX columns see [formats.md](formats.md).

---

## How the orchestrator runs

Both the CLI and the Streamlit Export tab call
`core.export.export_pipeline.run_export()` with three dataclasses:

- `ExportSelection` — what to export (parquets, comparison summary,
  mirror plots, summary plots).
- `PlotSettings` — the rendering knobs (`export_formats`,
  `export_resolution`, `figure_size`, `workers`,
  `plot_unpaired_queries`).
- `ExportCallbacks` — `on_step`, `on_plot_progress`, `on_warning` —
  hooks the front-end uses for progress UI.

Mirror plot tasks are built by
`core.pipeline.query_processing.build_plot_tasks` (one
`ComparisonTask` per query/reference pair) and dispatched through
`core.common.task_dispatch.run_tasks` to the
`core.plotting.plot_worker.process_single_task` worker. Summary
plots are rendered on the main process by
`core.plotting.summary_pipeline.generate_all_summary_plots`.

---

## Mirror plots

One plot per `(query spectrum, reference spectrum)` pair, rendered in
parallel by `plot_worker.process_single_task`. Each task carries
everything the worker needs to be picklable — the annotated query frame,
the reference frame (or `None` for single-spectrum tasks), labels,
metrics text, and per-format output paths — so workers do not need to
share any state with the main process.

### File formats

| Format | Directory | Renderer | Typical use |
|---|---|---|---|
| PNG | `png/` | Matplotlib `Agg` backend | Publication figures, slides |
| SVG | `svg/` | Matplotlib SVG backend | Vector editing in Illustrator / Inkscape |
| PDF | `pdf/` | Matplotlib PDF backend | Vector handoff with font embedding |
| HTML | `html/` | Plotly | Interactive zoom/hover/peak selection |

PNG/SVG/PDF are saved by `_save_task_figure` (per-format kwargs in
`_FORMAT_SAVE_KWARGS`); HTML is saved separately by `_save_task_html`
because Plotly renders directly without a Matplotlib figure. All four
share the same logical layout: query on top (positive intensity),
reference mirrored below (negative intensity), annotated peaks coloured
by ion series, metrics anchored in the lower-left corner.

### File naming

Plots are organized into per-peptide subdirectories (base sequence +
modification state). Within each subdirectory:

```
{query_peptide}_vs_{reference_peptide}_{ref_file_stem}_z{charge}.{ext}
```

For example:

```
png/ASDDFDQASY_unmodified/ASDDFDQASY_vs_ref_ProT001_z2.png
```

Filenames are passed through `core.export.export_directories.
sanitize_filename` which strips filesystem-unsafe characters and
truncates to 100 chars.

Single-spectrum plots (no reference available, Koina disabled, or
`plot_unpaired_queries: true`) use the simpler form:

```
png/ASDDFDQASY_unmodified/ASDDFDQASY_z2.png
```

### Per-plot contents

- **Axes:** m/z horizontal, intensity vertical. Intensity units match
  the `intensity_mode` setting — Relative = fraction of base peak, Raw
  = absolute counts.
- **Annotations:** ion labels (`y7`, `b3++`) drawn above each annotated
  peak. Colours read from `ION_COLORS` in
  `core/common/defaults.py` (mirrored in `settings.json` for
  Streamlit).
- **Title:** peptidoform string + charge state.
- **Metrics panel:** the 7 comparison metrics + context counts (peak
  totals, annotation counts, shared).

---

## Summary plots

Run-wide aggregations across every comparison in
`comparison_summary.csv`. Produced by
`core.plotting.summary_pipeline.generate_all_summary_plots` and
written to `{run}/summary/`.

The default metric is **R² binned** (`DEFAULTS.summary_metric`),
configured via the `summary_metric` YAML key for the CLI or the dropdown
on the Streamlit Export tab.

### Plot types

| File | Function | What it shows |
|---|---|---|
| `dual_axis.{png,svg,pdf}` | `summary_dual_axis.plot_classic_dual_axis` | Horizontal count bars (left axis, # comparisons per peptidoform) paired with metric boxplots (right axis, metric distribution) — one row per peptidoform |
| `ridgeline_r2_binned.{png,svg,pdf}` *or* `ridgeline.{png,svg,pdf}` | `summary_ridgeline.plot_ridgeline` | Per-peptidoform KDE ridgeline of the chosen metric. The `_r2_binned` suffix is used **only** when `metric == DEFAULTS.summary_metric`; any other metric drops the suffix |

Both plots share the same data-prep step (`prepare_summary_data`), which
loads the comparison CSV (or accepts an in-memory frame), constructs the
`peptidoform` label (`Modification + "/" + Charge`), and aggregates per
peptidoform: count, mean, median, std, q25/q75, 95% CI, dominant charge.

There is **no** "heatmap" summary plot — `_PLOT_REGISTRY` contains only
`dual_axis` and `ridgeline`. The openms_insight Heatmap component on the
Explore tab is interactive-only and is not saved to disk.

### Registry

Plot types are registered in `summary_pipeline._PLOT_REGISTRY` as a
straight `{name: callable}` mapping. To add a new plot type:

1. Add a new module `core/plotting/summary_<name>.py` exporting a
   `plot_<name>(raw_frame, aggregated_frame, metric, …) -> Figure` that
   matches the existing renderer signature.
2. Import it in `summary_pipeline.py` and add an entry to
   `_PLOT_REGISTRY`.
3. (Optional) Add the type's name to `DEFAULTS.summary_plot_types` so
   it's rendered by default.

The registry keeps renderer contracts consistent — every plot consumes
the same `(raw_frame, aggregated_frame, metric, export_resolution,
title)` keyword set.

### Customizing

- `metric` — any column listed in `METRIC_COLUMNS`
  (`core/plotting/_summary_shared.py`): `R² match`, `R² binned`,
  `Cosine`, `Entropy`, `Dice`, `Jaccard`, `SharedRatio`.
- `plot_types` — subset of `{dual_axis, ridgeline}` to render. Empty
  list → no summary plots.
- `formats` — list like `["png", "svg"]` (defaults to PNG only).
- `export_resolution` — DPI passed through to `Figure.savefig`
  (default 300).

Summary plots run on the main process; there's no parallel speed-up
because each call is O(total comparisons) and fast.

---

## Gallery & history

Two Streamlit-only post-export views let you browse past runs without
leaving the app. Both live in `streamlit_app/views/export/`.

### Plot gallery

Rendered by `gallery._render_plot_gallery(root_dir)`. After an export
completes, the Export tab shows an expandable tree:

- Outer expander: one per peptide-variant subdirectory.
- Inner expander: one per reference spectrum within that variant.
- Per-plot: inline image preview (PNG or SVG) plus a download button
  for each exported format.

The gallery prefers SVG over PNG when both exist (falls back to PNG if
SVG wasn't selected). HTML plots are listed but rendered as download
links, not embedded, to avoid Streamlit-iframe memory issues.

### Export history

The Export tab's history panel lists prior runs in the configured export
base directory, newest first. For each run you can:

- **View the fingerprint** (first 16 hex chars of the settings SHA-256).
- **Reload settings** — copies the run's `run_settings.json` back into
  the current workspace, so the Setup tab reflects that run's
  parameters.
- **Delete** — removes the numbered directory. Cannot be undone.

Matching runs (a run whose fingerprint equals the current Setup
parameters) are highlighted so you can tell at a glance whether a re-run
would reuse an existing directory or create a new one. See
[layout.md#re-run-behavior](layout.md#re-run-behavior) for the
fingerprint-matching rules.

---

## See also

- [layout.md](layout.md) — directory tree and `run_settings.json` schema.
- [formats.md](formats.md) — Parquet + CSV/XLSX column reference.
- [../metrics/comparison.md](../metrics/comparison.md) — what each
  metric in the summary CSV actually measures.
- [../architecture/core.md](../architecture/core.md) — the
  `core/plotting/` and `core/export/` subpackage tours.
