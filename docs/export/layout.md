# Export Directory Layout

What an export directory looks like on disk, how it's versioned, and what
`run_settings.json` contains. For file-level schemas (Parquet, CSV), see
[formats.md](formats.md). For plots, see [plots.md](plots.md).

Implementation references:

- Shared orchestrator (CLI + Streamlit): `mspecviz/export/runner.py`
- CLI thin wrapper: `cli/main.py:_run_export`
- Streamlit thin wrapper: `streamlit_app/views/export/pipeline.py`
- Versioned directory + fingerprint: `mspecviz/export/directories.py`

---

## Export directory structure

Each export run creates a versioned directory under the configured export
base path. The naming convention is:

```
{project_name}_{NNN}_{YYYYMMDD}/
```

For example: `my_project_001_20260327/`

- `{project_name}` is the sanitized project name from the Export tab or
  the `project_name` config key.
- `{NNN}` is a zero-padded sequence number, incrementing with each new
  configuration.
- `{YYYYMMDD}` is the date of the export run.

### Re-run behavior

MSpecViz fingerprints the run settings (a SHA-256 hash of all parameters
that affect output). When you re-run with identical settings:

- The existing directory is reused (renamed to today's date if needed).
- Plots already on disk are skipped, saving time.
- Only new or missing outputs are generated.
- `_last_run` in `run_settings.json` is bumped to the new timestamp.

When any setting changes, a new numbered directory is created.

### Directory layout

```
my_project_001_20260327/
    run_settings.json
    comparison_summary.csv
    comparison_summary.xlsx
    psm.parquet
    peaks.parquet
    seq.parquet
    heatmap.parquet
    summary/
        dual_axis_R2_binned.png
        ridgeline_R2_binned.png
        heatmap_all_metrics.png
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
for the formats you select. Within each format directory, plots are
organized into subdirectories by peptide variant (base sequence plus
modification state).

---

## `run_settings.json`

Each export directory contains a `run_settings.json` file that captures
all parameters used for the run. This enables reproducibility and
supports the export history features (reload settings, detect matching
runs).

| Field | Description |
|---|---|
| `_fingerprint` | SHA-256 hash (first 16 hex chars) of the canonical settings |
| `_created_at` | ISO timestamp of the initial export run |
| `_last_run` | ISO timestamp of the most recent re-run (updated on reuse) |
| `selected_idxmls` | List of input idXML file paths |
| `selected_mzmls` | List of input mzML file paths |
| `target_peptides` | Target peptide sequences |
| `reference_substring` | Reference file substring filter |
| `ion_types` | Fragment ion types used for annotation |
| `mz_tolerance` | m/z matching tolerance in Da |
| `add_neutral_losses` | Whether neutral losses were included |
| `add_first_prefix_ion` | Whether first prefix ions were included |
| `remove_precursor_peaks` | Whether precursor peaks were removed |
| `intensity_mode` | "Relative" or "Raw" |
| `koina_enabled` | Whether Koina AI prediction was active |
| `koina_model` | Koina model ID (if enabled) |
| `koina_collision_energy` | Collision energy (if enabled) |
| `koina_instrument_type` | Instrument type (if enabled) |
| `export_formats` | List of output formats (e.g., `["png", "svg"]`) |
| `export_resolution` | Raster resolution in DPI (dots per inch) |
| `figure_size` | Figure dimensions `[width, height]` |

The file is human-readable JSON — safe to edit by hand, though the
`_fingerprint` will no longer match after you do.
