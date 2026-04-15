# Data File Formats

Per-column schemas for every tabular artifact an export run produces.

- **Parquet** — four files, used both as the live in-app cache and as the
  portable data artifact. Built by `mspecviz.dataframes`.
- **CSV / XLSX** — one wide comparison summary per run, human-readable.

For plots see [plots.md](plots.md); for the directory layout see
[layout.md](layout.md).

---

## Parquet DataFrames

Four Parquet files are exported, corresponding to the four Polars
LazyFrames built for openms_insight. All use zstd compression.

### `psm.parquet` — PSM table

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

### `peaks.parquet` — Peak data

One row per peak in each loaded spectrum. Used by openms_insight LinePlot
and SequenceView.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | m/z value |
| `y` | Float32 | Intensity (normalized or raw, depending on intensity mode) |
| `highlight` | Int32 | 1 if the peak is annotated, 0 if unannotated |
| `annotation` | String | Ion label (e.g., "y7", "b3++") or "unannotated" |
| `ion_type` | String | Ion series (e.g., "b", "y", "a", "precursor", "unannotated") |
| `psm_id` | Int32 | Foreign key linking to the PSM table `id` column |

### `seq.parquet` — Sequence data

One row per PSM. Used by openms_insight SequenceView for rendering the
peptide sequence with modification markers.

| Column | Type | Description |
|---|---|---|
| `sequence` | String | Peptide in OpenMS parenthesis notation (e.g., "SHC(Carbamidomethyl)IAEVEK") |
| `precursor_charge` | Int32 | Precursor charge state |
| `psm_id` | Int32 | Foreign key for cross-component linking |

### `heatmap.parquet` — Heatmap data

One row per peak across all loaded spectra. Used by openms_insight
Heatmap for RT-by-m/z visualization.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | Retention time in seconds |
| `y` | Float32 | m/z value |
| `intensity` | Float32 | Peak intensity |

---

## CSV / XLSX comparison summary

Both `comparison_summary.csv` and `comparison_summary.xlsx` contain one
row per query-reference comparison. Columns are identical across the two
formats.

### Identification columns

| Column | Description |
|---|---|
| `query_file` | Source mzML filename of the query spectrum |
| `ref_file` | Source mzML filename of the reference spectrum (or `"Koina"` for an AI-predicted reference) |
| `peptide` | Peptidoform string |
| `modification` | Modification state |
| `charge` | Precursor charge |

### 7 comparison metrics

See [../metrics/comparison.md](../metrics/comparison.md) for definitions.

| Column | Description |
|---|---|
| `R2m` | R² matching |
| `R2b` | R² binned |
| `Cos` | Cosine similarity |
| `Ent` | Spectral entropy similarity |
| `Dice` | Dice coefficient |
| `Jacc` | Jaccard index |
| `ShR` | Shared peak ratio |

### Context counts

| Column | Description |
|---|---|
| `total_q` | Total peaks in query |
| `total_r` | Total peaks in reference |
| `ann_q` | Annotated peaks in query |
| `ann_r` | Annotated peaks in reference |
| `shared` | Shared annotated ions |

### 17 ion match statistics

All prefixed with `ims_`. See
[../metrics/ion-match-stats.md](../metrics/ion-match-stats.md) for full
descriptions.

| Column | Short description |
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

The XLSX variant auto-sizes columns and applies a freeze pane at row 1,
column B for readability. Otherwise it carries identical data.
