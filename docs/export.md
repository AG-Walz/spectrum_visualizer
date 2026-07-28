# Export

What an export run writes, how directories are versioned, and the file formats.
Both front ends call `core/export.py:run_export` with the same `ExportSelection` /
`PlotSettings` / `ExportCallbacks`, so CLI and Streamlit produce byte-identical
output. Per-column Parquet schemas live as self-documenting dicts in
`core/dataframes.py` — those constants are the source of truth, not this doc.

---

## Directory layout and versioning

Each run creates `{project_name}_{NNN}_{YYYYMMDD}/` under the configured base dir
(`project_name` sanitised and truncated to 100 chars; `NNN` a zero-padded sequence).

```
MSpecViz_analysis_001_20260327/
    run_settings.json
    comparison_summary.csv          # + .xlsx when openpyxl is installed
    psm.parquet  peaks.parquet  seq.parquet  heatmap.parquet
    summary/      ridgeline_r2_binned.png, multi_boxplot_r2_binned.png
    png/ svg/ pdf/ html/            # only the formats you selected
        <peptidoform>/              # one subdir per base sequence + modification
            <query>_vs_<ref>_z<charge>.<ext>
```

**Re-run behaviour.** A run is fingerprinted by a SHA-256 of its sorted
`run_settings` (first 16 hex chars → `_fingerprint`). On export,
`find_or_create_export_directory` scans the base dir: if a folder's `_fingerprint`
matches, it is renamed to today's date, reused, and its `_last_run` bumped —
already-rendered plots are skipped. No match → a new `NNN = max + 1` folder, created
only at `commit_export_directory` (so aborted runs leave nothing behind). To force a
fresh directory, change any field in `RUN_SETTINGS_FIELDS`.

**`run_settings.json`** captures the run for reproducibility and powers the
"reload settings" history feature. Three `_`-prefixed keys are its only home:
`_fingerprint`, `_created_at`, `_last_run`. The rest is exactly `RUN_SETTINGS_FIELDS`
(source of truth) plus `peptides` (added by `build_run_settings`). Editing it breaks
the fingerprint match, so the next run allocates a fresh directory.

---

## Data files

**Parquet** — the four LazyFrames openms_insight consumes, built by
`core/dataframes.py`: `psm.parquet` (one row per PSM, including its `retention_time`
in seconds — `0.0` when the identification carried none), `peaks.parquet` (one row per
peak, with ion annotations), `seq.parquet` (one row per PSM, OpenMS parenthesis
notation for SequenceView), `heatmap.parquet` (one row per peak: RT, m/z, intensity).
Schemas: the `*_SCHEMA` constants in that module.

**CSV / XLSX** — `comparison_summary.{csv,xlsx}`, one row per query-reference pair,
identical columns in both. The column set and order are owned by
`COMPARISON_SUMMARY_COLUMNS` in `core/pipeline/load_annotate.py` (built by
`build_comparison_frame`); read that constant for the exact list. It runs in three
groups: **identification** (`psm_id` — the FK to `psm.parquet` — plus the peptide,
charge, query/reference labels, spectrum IDs and the RT block `query_rt` / `ref_rt` /
`rt_delta`, the latter being query minus reference RT in seconds and `NaN` when the
reference is predicted and so has no retention time), the **7 metrics** under their
display names (`METRIC_KEY_TO_COLUMN_NAME` — see [metrics.md → the seven metrics](metrics.md#the-seven-metrics)), and
**context counts** (`query`/`reference` peak and annotation tallies plus `shared`).

The context counts are essential for sanity-checking a metric: a high cosine on
`shared = 2` is far weaker evidence than on `shared = 20` (see
[metrics.md → context columns](metrics.md#context-columns)).

---

## Plots

**Mirror plots** — one per `(query, reference)` pair, rendered in parallel by the
picklable `core/plotting/plot_worker.py:process_single_task`. Query on top (positive
intensity), reference mirrored below, peaks coloured by ion series, metrics in the
lower-left. Single-spectrum tasks (no reference) drop the `_vs_<ref>` part of the
name. Filenames pass through `sanitize_filename`.

Every plot spans `0 -> max(m/z)` of its **peptide group**, not of its own peaks:
`build_plot_tasks` resolves one bound per peptidoform via `peptidoform_mz_axis_maxima`
(over all of that peptide's queries and references) and stamps it on each task as
`mz_axis_max`. So the PNGs inside one peptide folder — and the same peptide's plot in
the Explore tab, which resolves the identical bound — can be flipped through without
the x-axis rescaling underneath (#114).

| Format | Renderer | Use |
|---|---|---|
| PNG | Matplotlib `Agg` | Publication figures, slides |
| SVG | Matplotlib SVG | Vector editing |
| PDF | Matplotlib PDF | Vector handoff with embedded fonts |
| HTML | Plotly | Interactive zoom/hover/peak selection |

**Summary plots** — run-wide aggregations in `{run}/summary/`, by
`core/plotting/summary.py:generate_all_summary_plots`. Only two types exist
(`SUMMARY_RENDERERS`): `ridgeline` (per-peptidoform KDE of the chosen metric) and
`multi_boxplot` (one horizontal boxplot per peptidoform with every data point). One
set is rendered per measure in `summary_metrics`, each saved as
`<plot_type>_<metric_slug>.<format>` (e.g. `multi_boxplot_cosine.png`).
There is no "heatmap" summary plot — the Explore-tab Heatmap is interactive-only and
never written to disk.

**Gallery & history (Streamlit-only)** — `views/export/`. The gallery shows an
expandable tree per peptidoform/reference with inline image previews (prefers SVG,
HTML as download links). The history panel lists prior runs newest-first; per run
you can view its fingerprint, reload its settings into the current workspace, or
delete it. Runs whose fingerprint matches the current Setup parameters are
highlighted (a re-run would reuse them).
