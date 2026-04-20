# Core Package: `mspecviz/`

Seven subpackages, each exposing a public API via a re-exporting
`__init__.py`. All cross-package imports should go through those
subpackage-level re-exports (`from mspecviz.search import ...`), never
reach into internal modules (`from mspecviz.search.workers import ...`).

Two small modules sit at the root of `mspecviz/` because they're shared
across multiple subpackages — covered at the [end of this page](#parent-level-helpers).

---

## `mspecviz/common/`

Pure helpers with no domain logic. Imported by everything.

| Module | What it holds |
|---|---|
| `constants.py` | `STATE_PREFIX`, the fallback `DEFAULTS` dict (last-resort values if neither `settings.json` nor workspace `params.json` supply a key), and plot tuning constants (`ANN_OFFSET`, `MZ_WINDOW_PADDING`) |
| `ion_utils.py` | `ION_COLORS` palette + `get_ion_color()` resolver (shared by both mirror renderers) |
| `normalization.py` | Intensity normalization helpers (relative-to-base-peak conversion) |
| `parallel.py` | `run_tasks()` — the single `ProcessPoolExecutor` orchestrator every worker file uses. Handles progress callbacks in the parent process so Streamlit progress bars stay attached |
| `text_utils.py` | `sanitize_filename`, `parse_peptide_list`, `format_time`, `get_nested` (dict-path accessor used by the config builder) |

---

## `mspecviz/search/`

idXML parsing and peptide search across one or more identification files.

| Module | What it holds |
|---|---|
| `models.py` | `PeptideMatch` (single PSM row), `VariantGroup` (base sequence + modification state + matches), `format_inline_sequence()` |
| `idxml_parsing.py` | Modification abbreviation tables, `extract_spectra_data_map()` (recovers idXML → mzML filename resolution), `mod_key_from_peptidoform()` |
| `workers.py` | Module-level, picklable worker functions: `_search_single_idxml`, `_scan_single_idxml_peptides`, `count_psms_per_peptide`, `_build_match_from_psm` |
| `orchestration.py` | Public entry points: `search_peptides_in_idxmls()`, `get_unique_peptides_from_idxmls()`, `summarize_search_results()` |

---

## `mspecviz/spectra/`

mzML streaming + Koina AI prediction. Both produce the same `LoadedSpectrum`
dataclass so the annotation layer doesn't need to know the origin.

| Module | What it holds |
|---|---|
| `loader_models.py` | `LoadedSpectrum` dataclass, `StreamingResult`, `LoadResult`, `spectrum_lookup_key()` |
| `mzml_native_ids.py` | Native-ID regex tags (`scan=`, `index=`, `spectrum=`) and `clean_spectrum_id()` normaliser. **Note:** `index=` IDs from timsTOF/Bruker are 1-based — the loader auto-detects and corrects |
| `mzml_lookup.py` | `open_indexed_mzml_fast()` + `find_spectrum_index()` (index-backed O(1) spectrum lookup) |
| `loader_workers.py` | Picklable `_load_spectra_from_mzml` worker |
| `loader_orchestration.py` | `load_spectra_streaming()` — the public API: streams spectra in parallel, matches them to PSMs, dispatches references/queries |
| `koina_prediction.py` | `predict_koina_spectrum()`, `KOINA_MODELS` catalogue, `get_koina_settings()`, `INSTRUMENT_TYPES` |

---

## `mspecviz/annotation/`

Theoretical fragment generation, peak annotation, and the 7 comparison
metrics.

| Module | What it holds |
|---|---|
| `peak_theory.py` | `peptidoform_to_aasequence()` (psm_utils → pyopenms), `generate_theoretical_fragments()`, `DEFAULT_ION_TYPES`, `ALL_ION_TYPES` |
| `peak_matching.py` | `annotate_spectrum()`, `prepare_spectrum_for_plot()` — the main entry point used by the pipeline — plus the Koina-specific variants `prepare_koina_spectrum_for_plot()` |
| `metrics_intensity.py` | Intensity-based metrics: `compute_r2_matching`, `compute_r2_binned`, `compute_cosine_similarity` |
| `metrics_set.py` | Set-based + entropy metrics: Dice, Jaccard, shared-peak ratio, spectral entropy; plus `compute_all_metrics()`, `format_metrics_text()`, and `METRIC_SHORT_TO_DISPLAY` |

See [../metrics/comparison.md](../metrics/comparison.md) for formulas,
ranges, and interpretation.

---

## `mspecviz/dataframes/`

Polars LazyFrames consumed by the openms_insight components.

| Module | What it holds |
|---|---|
| `psm_tables.py` | `PSM_TABLE_SCHEMA`, `build_psm_table()`, `validate_dataframe()`, `_enforce_schema()` |
| `visualization_frames.py` | `PEAKS_SCHEMA`, `SEQUENCE_SCHEMA`, `HEATMAP_SCHEMA`, `build_peaks_dataframe`, `build_sequence_dataframe`, `build_heatmap_dataframe`, `save_to_parquet`, `_proforma_to_openms_str` |

See [../export/formats.md](../export/formats.md) for the on-disk Parquet
schemas.

---

## `mspecviz/plotting/`

Both mirror renderers (Plotly and Matplotlib), three summary plot types,
and the parallel plot-generation worker.

| Module | What it holds |
|---|---|
| `_mirror_shared.py` | `_ensure_annotation_col`, `_prepare_spectrum_data`, `_get_annotated_mask` — shared between the two mirror renderers |
| `mirror_plotly.py` | `create_plotly_mirror_plot`, `create_plotly_single_plot`, `_add_plotly_spectrum_traces` — interactive HTML output |
| `mirror_matplotlib.py` | `create_mirror_plot`, `create_single_spectrum_plot`, `figure_to_bytes` — static PNG/SVG/PDF output |
| `summary_dual_axis.py` | `plot_classic_dual_axis` (PSM counts bar chart + metric box-plots) |
| `summary_ridgeline.py` | `plot_ridgeline` (KDE stacked by peptidoform) |
| `summary_orchestration.py` | `generate_all_summary_plots`, `prepare_summary_data`, `METRIC_COLUMNS`, `DEFAULT_METRIC`, `SUMMARY_COLORS`, `PLOT_TYPES`, `_PLOT_REGISTRY` |
| `worker_tasks.py` | `ComparisonTask`, `ComparisonResult` dataclasses (picklable work units) |
| `worker_execution.py` | `process_single_task`, `run_parallel_comparisons` — ProcessPoolExecutor entry point for batch plot export |

---

## `mspecviz/pipeline/`

The glue between search, load, annotate, compare, and plot. **This is the
layer both the CLI and Streamlit share** — any new feature that should run
identically from both front ends belongs here.

| Module | What it holds |
|---|---|
| `annotation_steps.py` | `_annotate_query_spectrum`, `_cache_reference_spectra`, `get_koina_prediction`, `_try_koina_fallback` |
| `orchestration.py` | `process_query_result` (called once per streamed spectrum), `build_comparison_row`, `build_plot_tasks`, `format_spectrum_labels`, `get_source_label`, `build_psm_to_variant_mapping` |

---

## Shared cross-cutting helpers

- `common/path_utils.py` — `disambiguate_paths()` keeps same-named files
  from different directories distinct in display labels. Used by both
  `search/` and `spectra/`.
- `export/directories.py` — versioned export directory management +
  settings fingerprinting (SHA-256 of canonical settings → stable 16-char
  hash). Consumed by `export/runner.py`, which is the single entry point
  for CLI + Streamlit exports.
