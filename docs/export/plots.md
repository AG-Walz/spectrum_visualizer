# Plot Outputs

What plot files an export run produces, how they are named, and the
Streamlit-only gallery/history views that let you browse past runs.

- For per-spectrum **mirror plots** see [Mirror plots](#mirror-plots).
- For per-run **summary plots** see [Summary plots](#summary-plots).
- For Streamlit's post-export browsers see [Gallery & history](#gallery--history).

For the directory tree see [layout.md](layout.md); for the
comparison-summary CSV/XLSX columns see [formats.md](formats.md).

---

## Mirror plots

One plot per `(query spectrum, reference spectrum)` pair, generated in
parallel by `mspecviz.plotting.worker_execution.run_parallel_comparisons`.
Each pair is a `ComparisonTask` that the worker pool renders into as many
formats as you selected.

### File formats

| Format | Directory | Renderer | Typical use |
|---|---|---|---|
| PNG | `png/` | Matplotlib `Agg` backend | Publication figures, slides |
| SVG | `svg/` | Matplotlib SVG backend | Vector editing in Illustrator / Inkscape |
| PDF | `pdf/` | Matplotlib PDF backend | Vector handoff with font embedding |
| HTML | `html/` | Plotly | Interactive zoom/hover/peak selection |

All four share the same layout: query on top (positive intensity), reference
mirrored below (negative intensity), annotated peaks coloured by ion series,
metrics panel anchored to the lower-left corner.

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

Disambiguation suffixes are appended when two references would produce the
same filename (see `mspecviz/path_utils.py`).

Single-spectrum plots (no reference available, Koina disabled) use the
simpler form:

```
png/ASDDFDQASY_unmodified/ASDDFDQASY_z2.png
```

### Per-plot contents

- **Axes:** m/z horizontal, intensity vertical. Intensity units match the
  `intensity_mode` setting — Relative = fraction of base peak, Raw = absolute
  counts.
- **Annotations:** ion labels (`y7`, `b3++`) drawn above each annotated
  peak. Colours read from `settings.json > ion_colors`.
- **Title:** peptidoform string + charge state.
- **Metrics panel:** the 7 comparison metrics + context counts
  (annotated-Q / annotated-R / shared).
- **Footer:** query file, reference file, tolerance, ion types used.

---

## Summary plots

Run-wide aggregations across every comparison in `comparison_summary.csv`.
Produced by `mspecviz.plotting.summary_orchestration.generate_all_summary_plots`
and written to `{run}/summary/`.

The default metric is **R² binned**, selected via the `metric` argument (CLI
flag `--summary-metric`, Streamlit dropdown in the Export tab).

### Plot types

| File | Function | What it shows |
|---|---|---|
| `dual_axis.{png,svg,pdf}` | `plot_classic_dual_axis` | Horizontal count bars (left axis, # comparisons per peptidoform) paired with metric boxplots (right axis, metric distribution) — one row per peptidoform |
| `ridgeline_r2_binned.{png,svg,pdf}` | `plot_ridgeline` | Per-peptidoform KDE ridgeline of the chosen metric. Visualizes distribution shape, including bimodality that boxplots hide |

Both plots share the same data-preparation step
(`prepare_summary_data`), which aggregates per-peptidoform counts, mean,
median, and quartiles of the selected metric.

### Registry

Plot types are registered in `summary_orchestration._PLOT_REGISTRY`. Adding
a new plot type means adding:

1. A new module under `mspecviz/plotting/summary_<name>.py` exporting a
   `plot_<name>(raw_df, agg_df, metric, …) -> Figure` function.
2. An entry in `_PLOT_REGISTRY` inside `generate_all_summary_plots`.
3. The string name into `PLOT_TYPES`.

The registry keeps renderer contracts consistent — every plot consumes the
same `(raw_df, agg_df, metric, figure_size, export_dpi, title)` signature.

### Customizing

- `metric` — one of `R² match`, `R² binned`, `Cosine`, `Entropy`, `Dice`,
  `Jaccard`, `SharedRatio` (see `METRIC_COLUMNS`).
- `plot_types` — subset of `PLOT_TYPES` to render (default: all).
- `formats` — list like `["png", "svg"]` (default: PNG only).
- `export_dpi` — passed through to Matplotlib (default: 300).

Summary plots are rendered on the main process only; there is no parallel
speed-up to gain because each call is O(total comparisons) and fast.

---

## Gallery & history

Two Streamlit-only post-export views let you browse past runs without
leaving the app. Both live in `streamlit_app/views/export/`.

### Plot gallery

Rendered by `gallery._render_plot_gallery(root_dir)`. After an export
completes, the Export tab shows an expandable tree:

- Outer expander: one per peptide-variant subdirectory.
- Inner expander: one per reference spectrum within that variant.
- Per-plot: inline image preview (PNG or SVG) plus a download button for
  each exported format.

The gallery prefers SVG over PNG when both exist (falls back to PNG if SVG
wasn't selected). HTML plots are listed but rendered as download links, not
embedded, to avoid Streamlit-iframe memory issues.

Only the first 20 plots per reference are shown inline — the remainder are
still on disk, just not auto-loaded.

### Export history

The Export tab's history panel lists prior runs in the configured export
base directory, newest first. For each run you can:

- **View the fingerprint** (first 16 hex chars of the settings SHA-256).
- **Reload settings** — copies the run's `run_settings.json` back into the
  current workspace, so the Setup tab reflects that run's parameters.
- **Delete** — removes the numbered directory. Cannot be undone.

Matching runs (a run whose fingerprint equals the current Setup parameters)
are highlighted so you can tell at a glance whether a re-run would reuse an
existing directory or create a new one. See
[layout.md#re-run-behavior](layout.md#re-run-behavior) for the
fingerprint-matching rules.

---

## See also

- [layout.md](layout.md) — directory tree and `run_settings.json` schema.
- [formats.md](formats.md) — Parquet + CSV/XLSX column reference.
- [../metrics/comparison.md](../metrics/comparison.md) — what each metric
  in the summary CSV actually measures.
- [../architecture/core.md#mspecvizplotting](../architecture/core.md#mspecvizplotting)
  — where the plot renderers live in the codebase.
