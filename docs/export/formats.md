# Data File Formats

Per-column schemas for every tabular artifact an export run produces.

- **Parquet** — four files, the same Polars LazyFrames the openms_insight
  Explorer consumes at runtime. Built by `core.dataframes`
  (`build_psm_table` + `build_visualization_frames`).
- **CSV / XLSX** — one wide comparison summary per run, human-readable.
  The frame is assembled by `core.pipeline.annotation_steps.
  build_comparison_row` and serialised by
  `core.export.export_pipeline._write_comparison_summary`.

For plots see [plots.md](plots.md); for the directory layout see
[layout.md](layout.md).

---

## Parquet DataFrames

Four Parquet files mirror the four LazyFrames built for openms_insight.

### `psm.parquet` — PSM table

One row per peptide-spectrum match. Schema constant: `PSM_TABLE_SCHEMA`
in `core/dataframes/psm_tables.py`.

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
and SequenceView. Schema: `PEAKS_SCHEMA` in
`core/dataframes/insight_frames.py`.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | m/z value |
| `y` | Float32 | Intensity (normalised or raw, depending on intensity mode) |
| `highlight` | Int32 | 1 if the peak is annotated, 0 if unannotated |
| `annotation` | String | Ion label (e.g., "y7", "b3++") or "unannotated" |
| `ion_type` | String | Ion series (e.g., "b", "y", "a", "precursor", "unannotated") |
| `psm_id` | Int32 | Foreign key linking to the PSM table `id` column |

### `seq.parquet` — Sequence data

One row per PSM. Used by openms_insight SequenceView for rendering the
peptide sequence with modification markers. Schema: `SEQUENCE_SCHEMA`.

| Column | Type | Description |
|---|---|---|
| `sequence` | String | Peptide in OpenMS parenthesis notation (e.g., "SHC(Carbamidomethyl)IAEVEK") |
| `precursor_charge` | Int32 | Precursor charge state |
| `psm_id` | Int32 | Foreign key for cross-component linking |

### `heatmap.parquet` — Heatmap data

One row per peak across all loaded spectra. Used by openms_insight
Heatmap for RT-by-m/z visualization. Schema: `HEATMAP_SCHEMA`.

| Column | Type | Description |
|---|---|---|
| `x` | Float32 | Retention time in seconds |
| `y` | Float32 | m/z value |
| `intensity` | Float32 | Peak intensity |

---

## CSV / XLSX comparison summary

Both `comparison_summary.csv` and `comparison_summary.xlsx` contain one
row per query-reference comparison. Columns are identical across the two
formats; the XLSX is a straight `pandas.to_excel(... engine="openpyxl")`
dump. If `openpyxl` is not installed the XLSX is skipped and a warning is
emitted.

### Identification columns

| Column | Description |
|---|---|
| `psm_id` | Foreign key matching `id` in `psm.parquet` |
| `Peptide` | Base sequence (uppercase, unmodified) |
| `Modification` | Modification state (empty string when variant_key == base_sequence) |
| `Charge` | Precursor charge |
| `Query` | Source label of the query spectrum (`<basename>:<spectrum_id>`) |
| `Reference` | Source label of the reference (or `Koina/<model>` for an AI-predicted reference) |
| `query_spectrum_id` | Native spectrum ID inside the query mzML |
| `ref_spectrum_id` | Native spectrum ID inside the reference mzML |
| `query_rt` | Query retention time (seconds) |
| `ref_rt` | Reference retention time (seconds) |

### 7 comparison metrics

Column names are the display strings from `METRIC_SHORT_TO_DISPLAY`
(`core/annotation/comparison_metrics.py`), not the short codes. See
[../metrics/comparison.md](../metrics/comparison.md) for definitions.

| Column | Short code | Description |
|---|---|---|
| `R² match` | R2m | R² (squared Pearson correlation) on shared annotated intensities |
| `R² binned` | R2b | R² on binned all-peak intensity profile |
| `Cosine` | Cos | Cosine similarity on shared annotated intensities |
| `Entropy` | Ent | Spectral entropy similarity (`ms_entropy`) |
| `Dice` | Dice | Dice coefficient on annotation sets |
| `Jaccard` | Jacc | Jaccard index on annotation sets |
| `SharedRatio` | ShR | Shared / smaller-set size |

### Context counts

| Column | Description |
|---|---|
| `query_peak_count` | Total peaks in the query spectrum |
| `reference_peak_count` | Total peaks in the reference spectrum |
| `query_annotation_count` | Number of **unique** annotation labels in the query (duplicates collapse) |
| `reference_annotation_count` | Number of **unique** annotation labels in the reference |
| `shared` | Number of annotation labels present in both spectra (\|S\|) |

These counts are essential for sanity-checking metric values — see
[../metrics/interpretation.md](../metrics/interpretation.md#cross-checks).
A high cosine on `shared = 2` is much weaker evidence than on
`shared = 20`.
