# Core Package: `core/`

Eight subpackages, each exposing a public API via a re-exporting
`__init__.py`. Cross-package imports should always go through those
subpackage-level re-exports (`from core.search import ...`), never
reach into internal modules (`from core.search.idxml_workers import ...`).

---

## Top-level modules

Pure helpers with no domain logic. Imported by everything.

| Module | What it holds |
|---|---|
| `defaults.py` | The frozen `Defaults` dataclass + `DEFAULTS` instance (single source of truth for every tunable parameter), the `KOINA_MODELS` catalogue, `INSTRUMENT_TYPES`, `ION_COLORS` (Material-Design palette), and `IMAGE_FORMATS` |
| `helpers.py` | `parse_peptide_list`, `format_time`, `ensure_unique_basenames` (raises `DuplicateBasenameError` when two selected files share a basename across directories), and `run_tasks()` — the single dispatcher every parallel worker uses. Wraps `joblib.Parallel` and yields `(task, ok, result)` triples so callers can update progress bars and surface per-file errors in the parent process. Falls back to a sequential loop when `workers=1` |

`STATE_PREFIX` (the Streamlit session-state namespace) lives in
`streamlit_app/utils/state/state_core.py`, not here — it's a UI-layer concern.

---

## `core/search/`

idXML parsing and peptide search across one or more identification files.

| Module | What it holds |
|---|---|
| `idxml_workers.py` | Module-level, picklable workers: `search_single_idxml` (returns matches + per-file errors), `count_psms_per_peptide`. Both are dispatched by `core.helpers.run_tasks` |
| `peptide_search.py` | Public entry points: `search_peptides_in_idxmls()`, `count_psms_in_idxmls()`, `summarize_search_results()`. Wires the optional `progress_callback` / `error_callback` into the dispatcher so the CLI's tqdm bar and Streamlit's progress widget can attach |

Public API: `from core.search import search_peptides_in_idxmls,
count_psms_in_idxmls, summarize_search_results`.

---

## `core/spectra/`

mzML streaming + Koina AI prediction. Both produce the same `LoadedSpectrum`
dataclass so the annotation layer doesn't need to know the spectrum's origin.

| Module | What it holds |
|---|---|
| `loaded_spectrum.py` | `LoadedSpectrum` dataclass + `StreamingResult` / `LoadResult` containers, `spectrum_lookup_key()`, `source_label()` |
| `mzml_native_ids.py` | `extract_id_number()` and `clean_spectrum_id()` for normalising idXML/mzML native IDs (`scan=`, `index=`, `spectrum=`). **timsTOF/Bruker `index=` IDs are 1-based** — the loader auto-detects and shifts |
| `mzml_lookup.py` | `open_indexed_mzml_fast()` + `find_spectrum_index()` — index-backed O(1) spectrum lookup, with `_detect_id_format` working out the file's ID scheme on first access |
| `mzml_workers.py` | Picklable `_load_spectra_from_mzml` worker plus precursor-resolution helpers |
| `spectrum_streaming.py` | `load_spectra_streaming()` — the public API. Builds an mzML → matches map, runs reference loads first (so queries can be paired during streaming), then yields query results lazily |
| `koina_prediction.py` | `predict_koina_spectrum()`, the 14-entry `KOINA_MODELS` catalogue, `get_koina_settings()`, `INSTRUMENT_TYPES`, and `convert_peptide_to_unimod()` |

---

## `core/annotation/`

Theoretical fragment generation, peak annotation, and the 7 comparison
metrics. Theoretical fragments are produced by `spectrum_utils.proforma`
inside `annotate_spectrum()` — there is no separate `peak_theory` module.

| Module | What it holds |
|---|---|
| `peak_matching.py` | `annotate_spectrum()` (proforma → annotated `MsmsSpectrum`), `prepare_spectrum_for_plot()` (main entry point used by the pipeline), `prepare_koina_spectrum_for_plot()` (Koina-output variant), `get_ion_type()` for bucket resolution |
| `intensity_metrics.py` | Intensity-based metrics: `compute_r2_matching` (R²m), `compute_r2_binned` (R²b), `compute_cosine_similarity` (Cos). All take aligned query/reference DataFrames; return NaN on zero-norm or <2 shared points |
| `comparison_metrics.py` | `compute_all_metrics()` (the public batch entry point — returns R²m, R²b, Cos, Ent, Dice, Jacc, ShR + 5 context counts), `compute_spectral_entropy_similarity()` (wraps `ms_entropy` with a 0.02-Da tolerance), `format_metrics_text()`, `METRIC_SHORT_TO_DISPLAY` |

See [../metrics/comparison.md](../metrics/comparison.md) for formulas,
ranges, and interpretation.

---

## `core/dataframes/`

Polars LazyFrames consumed by the openms_insight components.

| Module | What it holds |
|---|---|
| `psm_tables.py` | `PSM_TABLE_SCHEMA`, `build_psm_table()` — flattens variant groups into one row per PSM and returns `(frame, psm_lookup)` |
| `insight_frames.py` | `PEAKS_SCHEMA`, `SEQUENCE_SCHEMA`, `HEATMAP_SCHEMA` plus the builders `build_peaks_dataframe`, `build_sequence_dataframe`, `build_heatmap_dataframe`, and the convenience `build_visualization_frames` that returns all three together |

See [../export/formats.md](../export/formats.md) for the on-disk Parquet
schemas.

---

## `core/plotting/`

A single mirror renderer (Plotly + Matplotlib via a shared backend), two
summary plot types (dual-axis + ridgeline), and the parallel plot-generation
worker.

| Module | What it holds |
|---|---|
| `mirror.py` | `create_spectrum_plot()` — hand-rolled renderer that branches on `backend="ms_matplotlib"` vs `"ms_plotly"`. Pyopenms-viz' `plot_spectrum` was abandoned because it scrambles ion-type colouring when ion types are interleaved. Internal helpers `_draw_matplotlib` / `_draw_plotly` keep the rendering paths separate but share data prep |
| `plot_task_dataclasses.py` | `ComparisonTask`, `ComparisonResult` dataclasses — picklable work units kept separate from the worker module so they stay lightweight to import in main-process code |
| `plot_worker.py` | `process_single_task()` — runs in a worker process: annotates the query/reference pair, renders the plot, writes raster/vector/HTML outputs. Format-specific save kwargs live in `_FORMAT_SAVE_KWARGS` |
| `summary_dual_axis.py` | `plot_classic_dual_axis()` — PSM-counts bar chart + metric box-plot |
| `summary_ridgeline.py` | `plot_ridgeline()` — KDE per peptidoform, stacked |
| `summary_pipeline.py` | `generate_all_summary_plots()`, `prepare_summary_data()`, and the `_PLOT_REGISTRY` mapping plot-type names to the rendering functions above |
| `_summary_shared.py` | `METRIC_COLUMNS`, `SUMMARY_COLORS`, `_METRIC_LABELS`, plus shared helpers (`_auto_height`, `_style_ax`, `_subsample_group`) used by both summary renderers |

A "heatmap" summary plot type is **not** generated — only `dual_axis` and
`ridgeline`. The openms_insight Heatmap component on the Explore tab is a
separate, interactive widget driven by `HEATMAP_SCHEMA`.

---

## `core/pipeline/`

The glue between search, load, annotate, compare, and plot. **This is the
layer both the CLI and Streamlit share** — any new feature that should run
identically from both front ends belongs here.

| Module | What it holds |
|---|---|
| `kwargs_builders.py` | `build_annotation_kwargs()` and `build_koina_kwargs()` — take a generic `getter(key, default)` callable so the CLI (`SimpleNamespace`) and Streamlit (`params` dict) can both supply parameters without the core knowing the storage format |
| `annotation_steps.py` | The per-result helpers: `_annotate_query_spectrum`, `_cache_reference_spectra`, `_try_koina_fallback`, `get_koina_prediction`, `build_comparison_row` |
| `query_processing.py` | `process_query_result()` (one call per streamed spectrum, returns annotated query + metrics + plot-task seed), `build_plot_tasks()` (post-pass: turns annotated queries into `ComparisonTask` list), `build_psm_to_variant_mapping()`, `format_spectrum_labels()` |
| `streaming_orchestration.py` | `run_streaming_annotation()` — the top-level loop that drives `load_spectra_streaming` and applies `process_query_result` to each yielded result. `make_query_ready_callback()` builds the per-query callback the streamer calls in the main process |

---

## `core/export/`

Versioned export directories + shared export orchestration. Both the CLI and
the Streamlit Export tab call `run_export()` with their own progress
callbacks; every file-layout decision lives here so the two frontends produce
byte-identical outputs from the same inputs.

| Module | What it holds |
|---|---|
| `export_directories.py` | `find_or_create_export_directory()` + `commit_export_directory()` — the versioning + SHA-256 settings-fingerprint logic. Reuses an existing directory when the fingerprint matches; appends `_v2`, `_v3`, … otherwise. `sanitize_filename()` lives here too |
| `export_pipeline.py` | `run_export()` (entry point), the dataclasses passed in (`ExportSelection`, `PlotSettings`, `ExportCallbacks`, `ExportResult`), plus the orchestrated steps: parquet writes, comparison CSV/XLSX, parallel plot dispatch, summary plot generation |
| `run_settings.py` | `RUN_SETTINGS_FIELDS` (the single source of truth for which parameters are serialised) and `build_run_settings()` — same `getter` pattern as `kwargs_builders`, so the CLI and Streamlit hand in their own storage adapter |

See [../export/layout.md](../export/layout.md) for the directory layout and
[../export/formats.md](../export/formats.md) for the file schemas.
