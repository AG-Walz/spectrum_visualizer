# Export Directory Layout

What an export directory looks like on disk, how it's versioned, and what
`run_settings.json` contains. For file-level schemas (Parquet, CSV), see
[formats.md](formats.md). For plots, see [plots.md](plots.md).

Implementation references:

- Public orchestrator (CLI + Streamlit both call it):
  `core/export/export_pipeline.py::run_export()`
- Versioned directory + fingerprint logic:
  `core/export/export_directories.py` (`find_or_create_export_directory`,
  `commit_export_directory`, `sanitize_filename`)
- `run_settings.json` field list:
  `core/export/run_settings.py` (`RUN_SETTINGS_FIELDS`,
  `build_run_settings`)

The CLI wires it up in `cli/main.py::_run_export`; the Streamlit Export
tab does the same in `streamlit_app/views/export/export_actions.py`.
Both pass an `ExportSelection`, `PlotSettings`, and `ExportCallbacks`
dataclass — see [plots.md](plots.md#how-the-orchestrator-runs) for the
sequence.

---

## Export directory structure

Each export run creates a versioned directory under the configured export
base path. The naming convention is:

```
{project_name}_{NNN}_{YYYYMMDD}/
```

For example: `MSpecViz_analysis_001_20260327/`

- `{project_name}` is sanitised by `sanitize_filename` (strips
  filesystem-unsafe characters and truncates to 100 chars).
- `{NNN}` is a zero-padded sequence number, incrementing whenever a new
  fingerprint is seen under `{base_dir}`.
- `{YYYYMMDD}` is the date — see "Re-run behavior" below for how it
  updates.

### Re-run behavior

MSpecViz fingerprints the run settings: `_settings_fingerprint` is a
SHA-256 hash of `json.dumps(run_settings, sort_keys=True)` truncated to
the first 16 hex characters. This 16-char prefix is what lives in
`_fingerprint` inside `run_settings.json`.

On every export, `find_or_create_export_directory`:

1. Walks `{base_dir}` for folders matching `{project}_{NNN}_{YYYYMMDD}`.
2. Tracks the highest `NNN` it sees.
3. Reads each `run_settings.json` to find one whose `_fingerprint`
   matches the current settings.

If a match is found, the existing folder is renamed to today's date
(keeping its `NNN`) and reused — its `_last_run` timestamp is bumped on
write. Plots already on disk are skipped by the per-task `skip_plot`
check in the plot worker, so re-runs only generate missing outputs.

If no match is found, a new folder with `NNN = max + 1` is allocated.
Nothing is actually created until `commit_export_directory` runs (so
aborted runs don't leave empty folders).

To force a fresh directory, change any setting that participates in the
fingerprint (any field listed in `RUN_SETTINGS_FIELDS` — see below).

### Directory layout

```
MSpecViz_analysis_001_20260327/
    run_settings.json
    comparison_summary.csv
    comparison_summary.xlsx           # only when openpyxl is installed
    psm.parquet
    peaks.parquet
    seq.parquet
    heatmap.parquet
    summary/
        dual_axis.png                 # or .svg/.pdf depending on export_formats
        ridgeline_r2_binned.png       # ridgeline + default metric uses this name;
                                      # other metrics produce ridgeline.{ext}
    png/
        ASDDFDQASY_unmodified/
            ASDDFDQASY_vs_ref_ProT001_z2.png
            ...
        RTEFNLNQY_unmodified/
            ...
    svg/
        ASDDFDQASY_unmodified/
            ...
    pdf/
        ASDDFDQASY_unmodified/
            ...
    html/
        ASDDFDQASY_unmodified/
            ASDDFDQASY_vs_ref_ProT001_z2.html
            ...
```

Format subdirectories (`png/`, `svg/`, `pdf/`, `html/`) are created only
for the formats you select via `export_formats`. Within each format
directory, plots are organised into subdirectories by peptide variant
(base sequence plus modification state).

Note: there is **no** `heatmap_all_metrics.png`. Only two summary plot
types exist (`dual_axis` and `ridgeline`); the openms_insight Heatmap
component on the Explore tab is a separate, interactive widget that does
not produce a static export.

---

## `run_settings.json`

Each export directory contains a `run_settings.json` capturing the
parameters used for the run. It enables reproducibility and powers the
"reload settings" feature on the Streamlit Export-history panel.

The serialised fields are the single source of truth in
`core/export/run_settings.py::RUN_SETTINGS_FIELDS`. Both the CLI and
the Streamlit Export tab hand `build_run_settings()` a `getter(key,
default)` callable that knows how to read the value from their local
config representation; the field list itself is identical.

| Field | Description |
|---|---|
| `_fingerprint` | SHA-256 hash (first 16 hex chars) of the canonical settings dict |
| `_created_at` | ISO timestamp of the initial export run |
| `_last_run` | ISO timestamp of the most recent re-run (bumped on reuse) |
| `idxml_files` | List of input idXML file paths (post-glob expansion) |
| `mzml_files` | List of input mzML file paths (post-glob expansion) |
| `peptides` | Target peptide sequences (resolved list, not the YAML `peptide_source`) |
| `reference_substring` | Reference file substring filter |
| `ion_types` | Fragment ion types used for annotation |
| `mz_tolerance` | m/z matching tolerance in Da |
| `add_neutral_losses` | Whether H₂O / NH₃ neutral losses were included |
| `intensity_mode` | `"Relative"` or `"Raw"` |
| `koina_enabled` | Whether Koina AI prediction was active |
| `koina_model` | Koina model ID |
| `koina_collision_energy` | Collision energy (eV) |
| `koina_instrument_type` | Instrument type (`QE` / `LUMOS` / `TIMSTOF` / `SCIEXTOF`) |
| `export_formats` | List of output formats (`["png", "svg", "pdf", "html"]` subset) |
| `export_resolution` | Raster resolution in DPI |
| `figure_size` | Figure dimensions `[width, height]` in inches |
| `summary_metric` | Metric the summary plots are coloured/sorted by |
| `summary_plot_types` | Subset of `{dual_axis, ridgeline}` rendered |

The file is human-readable JSON — safe to inspect by hand. Editing it
will break the `_fingerprint` match, so a subsequent re-run will allocate
a fresh directory rather than reuse this one.
