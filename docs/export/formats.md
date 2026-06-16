# Data File Formats

Per-column schemas for every tabular artifact an export run produces.

- **Parquet** — four files, the same Polars LazyFrames the openms_insight
  Explorer consumes at runtime. Built by `core/dataframes.py`
  (`build_psm_table` + `build_visualization_frames`).
- **CSV / XLSX** — one wide comparison summary per run, human-readable.
  The frame's rows are assembled inline in
  `core/pipeline/load_annotate.py` and serialised inline in
  `core/export.py::run_export` (`to_csv` / `to_excel`).

For plots see [plots.md](plots.md); for the directory layout see
[layout.md](layout.md).

---

## Parquet DataFrames

Four Parquet files mirror the four LazyFrames built for openms_insight.
Each column schema lives as a self-documenting dict (with inline comments)
in `core/dataframes.py` — that constant is the source of truth, not this
doc.

- **`psm.parquet`** — one row per peptide-spectrum match (sequence,
  modification, charge, score, source files). Schema: `PSM_TABLE_SCHEMA`.
- **`peaks.parquet`** — one row per peak in each loaded spectrum, with ion
  annotations; consumed by openms_insight LinePlot and SequenceView.
  Schema: `PEAKS_SCHEMA`.
- **`seq.parquet`** — one row per PSM, the peptide in OpenMS parenthesis
  notation for SequenceView rendering. Schema: `SEQUENCE_SCHEMA`.
- **`heatmap.parquet`** — one row per peak across all spectra (RT, m/z,
  intensity) for the openms_insight Heatmap. Schema: `HEATMAP_SCHEMA`.

---

## CSV / XLSX comparison summary

Both `comparison_summary.csv` and `comparison_summary.xlsx` contain one
row per query-reference comparison. Columns are identical across the two
formats; the XLSX is a straight `pandas.to_excel(... engine="openpyxl")`
dump. If `openpyxl` is not installed the XLSX is skipped and a warning is
emitted.

These columns are assembled inline in
`core/pipeline/load_annotate.py` (plus the metric columns from
`core/comparison_metrics.py`) with no named constant — this doc is their
only consolidated reference.

### Identification columns

| Column | Description |
|---|---|
| `psm_id` | Foreign key matching `id` in `psm.parquet` |
| `Peptide` | Base sequence (uppercase, unmodified) |
| `Modification` | Modification state (empty string when peptidoform == base_sequence) |
| `Charge` | Precursor charge |
| `Query` | Source label of the query spectrum (`<basename>:<spectrum_id>`) |
| `Reference` | Source label of the reference (or `Koina/<model>` for an AI-predicted reference) |
| `query_spectrum_id` | Native spectrum ID inside the query mzML |
| `ref_spectrum_id` | Native spectrum ID inside the reference mzML |
| `query_rt` | Query retention time (seconds) |
| `ref_rt` | Reference retention time (seconds) |

### 7 comparison metrics

The seven metric columns use the display strings (not the short codes);
the short→display map is `METRIC_SHORT_TO_DISPLAY` in
`core/comparison_metrics.py`. See
[../metrics/comparison.md](../metrics/comparison.md) for what each metric
measures.

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
