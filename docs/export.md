# Export Reference

This document describes the export directory structure, file formats, and
data schemas produced by MSpecViz's Export tab.

---

## Export Directory Structure

Each export run creates a versioned directory under the configured export
base path. The naming convention is:

```
{project_name}_{NNN}_{YYYYMMDD}/
```

For example: `my_project_001_20260327/`

- `{project_name}` is the sanitized project name from the Export tab.
- `{NNN}` is a zero-padded sequence number, incrementing with each new
  configuration.
- `{YYYYMMDD}` is the date of the export run.

### Re-run Behavior

MSpecViz fingerprints the run settings (a SHA-256 hash of all parameters
that affect output). When you re-run with identical settings:

- The existing directory is reused (renamed to today's date if needed).
- Plots already on disk are skipped, saving time.
- Only new or missing outputs are generated.

When any setting changes, a new numbered directory is created.

### Directory Layout

```
my_project_001_20260327/
    run_settings.json
    comparison_summary.csv
    comparison_summary.xlsx
    psm.parquet
    peaks.parquet
    seq.parquet
    heatmap.parquet
    summary_plots/
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

Format subdirectories (png/, svg/, pdf/) are created only for the formats
you select. Within each format directory, plots are organized into
subdirectories by peptide variant (base sequence plus modification state).

---

## run_settings.json

Each export directory contains a `run_settings.json` file that captures all
parameters used for the run. This enables reproducibility and supports the
export history features (reload settings, detect matching runs).

Contents include:

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
| `add_losses` | Whether neutral losses were included |
| `add_first_prefix_ion` | Whether first prefix ions were included |
| `remove_precursor_peaks` | Whether precursor peaks were removed |
| `intensity_mode` | "Relative" or "Raw" |
| `koina_enabled` | Whether Koina AI prediction was active |
| `koina_model` | Koina model ID (if enabled) |
| `koina_ce` | Collision energy (if enabled) |
| `koina_instrument` | Instrument type (if enabled) |
| `export_formats` | List of output formats (e.g., ["png", "svg"]) |
| `dpi` | Image resolution |
| `figsize` | Figure dimensions [width, height] |

---

## Parquet DataFrames

Four Parquet files are exported, corresponding to the four Polars LazyFrames
built for openms_insight. All use zstd compression.

### psm.parquet -- PSM Table

One row per peptide-spectrum match.

| Column | Type | Description |
|---|---|---|
| `id` | Int32 | Unique PSM identifier (sequential, starting at 0) |
| `variant_key` | String | Peptide variant group key |
| `base_sequence` | String | Unmodified peptide sequence |
| `modification` | String | Modification state (e.g., "Oxidation" or "(unmodified)") |
| `display_label` | String | Human-readable label for the variant |
| `peptidoform` | String | Full ProForma peptidoform string |
| `charge` | Int32 | Precursor charge state |
| `score` | Float32 | Search engine score |
| `source_file` | String | mzML source filename |
| `spectrum_id` | String | Spectrum identifier within the mzML file |
| `idxml_file` | String | Source idXML filename |

### peaks.parquet -- Peak Data

One row per peak in each loaded spectrum. Used by openms_insight LinePlot and
SequenceView.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | m/z value |
| `y` | Float32 | Intensity (normalized or raw, depending on intensity mode) |
| `highlight` | Int32 | 1 if the peak is annotated, 0 if unannotated |
| `annotation` | String | Ion label (e.g., "y7", "b3++") or "unannotated" |
| `ion_type` | String | Ion series (e.g., "b", "y", "a", "precursor", "unannotated") |
| `psm_id` | Int32 | Foreign key linking to the PSM table `id` column |

### seq.parquet -- Sequence Data

One row per PSM. Used by openms_insight SequenceView for rendering the peptide
sequence with modification markers.

| Column | Type | Description |
|---|---|---|
| `sequence` | String | Peptide in OpenMS parenthesis notation (e.g., "SHC(Carbamidomethyl)IAEVEK") |
| `precursor_charge` | Int32 | Precursor charge state |
| `psm_id` | Int32 | Foreign key for cross-component linking |

### heatmap.parquet -- Heatmap Data

One row per peak across all loaded spectra. Used by openms_insight Heatmap for
RT-by-m/z visualization.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | Retention time in seconds |
| `y` | Float32 | m/z value |
| `intensity` | Float32 | Peak intensity |

---

## CSV/XLSX Comparison Summary

The comparison summary files (`comparison_summary.csv` and
`comparison_summary.xlsx`) contain one row per query-reference comparison
with the following columns:

### Identification Columns

| Column | Description |
|---|---|
| `query_file` | Source mzML filename of the query spectrum |
| `ref_file` | Source mzML filename of the reference spectrum (or "Koina") |
| `peptide` | Peptidoform string |
| `modification` | Modification state |
| `charge` | Precursor charge |

### 7 Comparison Metrics

| Column | Description |
|---|---|
| `R2m` | R-squared matching |
| `R2b` | R-squared binned |
| `Cos` | Cosine similarity |
| `Ent` | Spectral entropy similarity |
| `Dice` | Dice coefficient |
| `Jacc` | Jaccard index |
| `ShR` | Shared peak ratio |

### Context Counts

| Column | Description |
|---|---|
| `total_q` | Total peaks in query |
| `total_r` | Total peaks in reference |
| `ann_q` | Annotated peaks in query |
| `ann_r` | Annotated peaks in reference |
| `shared` | Shared annotated ions |

### 17 Ion Match Statistics

All prefixed with `ims_`. See [METRICS.md](METRICS.md) for full descriptions.

| Column | Description |
|---|---|
| `ims_NTermIonCurrentRatio` | N-terminal ion current fraction |
| `ims_CTermIonCurrentRatio` | C-terminal ion current fraction |
| `ims_TIC` | Total ion current |
| `ims_MIC` | Max ion current |
| `ims_annotated_tic_fraction` | Annotated TIC fraction |
| `ims_matched_prefix_ions_fraction` | Matched N-terminal ions fraction |
| `ims_matched_suffix_ions_fraction` | Matched C-terminal ions fraction |
| `ims_max_series_type` | Longest series ion type |
| `ims_max_series_charge` | Longest series charge |
| `ims_max_series_size` | Longest series length |
| `ims_coverage` | Backbone cleavage coverage |
| `ims_annotated_ion_fraction` | Annotated ion fraction |
| `ims_ms2_scan_index` | MS2 scan index |
| `ims_precursor_mz_error_ppm` | Precursor m/z error (ppm) |
| `ims_fragment_mass_error_mean_ppm` | Mean fragment mass error (ppm) |
| `ims_fragment_mass_error_stdev_ppm` | Fragment mass error std dev (ppm) |
| `ims_median_fragment_error_ppm` | Median fragment error (ppm) |

---

## Mirror Plots

### Static Plots (PNG, SVG, PDF)

Matplotlib-generated publication-quality mirror plots. Each plot shows:

- Upper panel: query spectrum with peaks pointing upward.
- Lower panel: reference spectrum with peaks pointing downward.
- Color-coded peaks by ion type (configurable via `ion_colors` in settings.json).
- Annotation labels on matched peaks.
- Metrics overlay text showing all 7 comparison metrics.
- Title with peptide, modification, charge, and source file information.

File naming: `{peptide}_{vs_ref_or_single}_{source}_{charge}.{format}`

Long filenames are truncated to 30 characters with a SHA-256 hash suffix to
guarantee uniqueness.

Default settings: 300 DPI, 12x8 inch figure size (configurable in the
Export tab or settings.json).

### Interactive HTML Plots

Plotly-generated interactive plots saved as standalone HTML files. These use
`include_plotlyjs="cdn"` to keep file sizes small (approximately 10 KB each
rather than 3 MB with embedded Plotly.js).

Interactive features include:

- Hover tooltips showing m/z, intensity, and annotation for each peak.
- Zoom, pan, and selection tools.
- Downloadable as PNG from the Plotly toolbar.

---

## Summary Analysis Plots

Three types of summary visualizations are generated from the comparison data,
saved in the `summary_plots/` subdirectory:

### Dual-Axis Plot

Bar chart of PSM counts per peptidoform overlaid with box-and-whisker plots
of a selected metric (default: R-squared binned). Individual data points are
scattered on top when group sizes are manageable (up to 200 per group).

### Ridgeline Plot

Kernel density estimate (KDE) distributions of a selected metric, stacked
vertically by peptidoform. A count sidebar shows the number of comparisons
per group.

### Multi-Metric Heatmap

Color-coded grid showing median values of all 7 metrics for each peptidoform.
Uses a red-yellow-green color scale where green indicates higher similarity.

---

## Export History

The Export tab maintains a history of all export runs found in the export
base directory. For each past run you can:

- **View settings** -- Inspect the run_settings.json to see what parameters
  were used.
- **Reload settings** -- Copy settings from a past run into the current Setup
  tab configuration.
- **Delete** -- Remove the export directory and all its contents.
- **ZIP download** -- Download the entire export directory as a ZIP archive
  directly from the browser.

---

## Plot Gallery

The plot gallery in the Export tab provides a browser for exported plots
organized by peptide variant and reference. Individual plots can be viewed
inline and downloaded without navigating the filesystem. The gallery reads
from the most recent export directory and supports all exported image formats.
