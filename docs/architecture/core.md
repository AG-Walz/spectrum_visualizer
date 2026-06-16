# Core Package: `core/`

Pure Python, no UI imports. Flat top-level modules plus three subpackages
(`pipeline/`, `plotting/`, `spectra/`). Each subpackage re-exports its public
API via `__init__.py`; import through those re-exports
(`from core.spectra import ...`), never reach in for private members.

---

## Top-level modules

| Module | What it holds |
|---|---|
| `defaults.py` | Frozen `Defaults` dataclass + `DEFAULTS` instance (single source of truth for every tunable), the `KOINA_MODELS` catalogue, `INSTRUMENT_TYPES`, and `ION_COLORS` (Material-Design palette) |
| `helpers.py` | `parse_peptide_list`, `format_time`, `ensure_unique_basenames` (raises `DuplicateBasenameError` when two selected files share a basename), and `run_tasks()` — the single parallel dispatcher. Wraps `joblib.Parallel` and yields `(task, ok, result)` triples so callers update progress bars and surface per-file errors in the parent process; sequential loop when `workers=1` |
| `peptide_lookup.py` | The only place idXML is parsed. `load_idxmls()` reads each file once into an in-memory `PSMList` (via the picklable worker dispatched through `run_tasks`), then `lookup_peptides()` and `count_psms_per_peptide()` scan that loaded list — the same file is never read twice |
| `peak_annotation.py` | Theoretical fragments via `spectrum_utils.proforma` (no separate peak-theory module). See module docstring for the per-symbol entry points (`prepare_spectrum_for_plot`, the Koina variant, `get_ion_type`) |
| `comparison_metrics.py` | `compute_all_metrics()` — the 7 metrics (R²m, R²b, Cos, Ent, Dice, Jacc, ShR) plus 5 context counts; all return NaN on zero-norm or <2 shared points. `METRIC_SHORT_TO_DISPLAY` maps short codes to labels. See [../metrics/comparison.md](../metrics/comparison.md) for formulas |
| `dataframes.py` | Polars-LazyFrame builders consumed by openms_insight: `build_psm_table` (one row per PSM, returns `(frame, psm_lookup)`) and `build_visualization_frames` (peaks + sequence + heatmap together), with the four schemas `PSM_TABLE_SCHEMA`, `PEAKS_SCHEMA`, `SEQUENCE_SCHEMA`, `HEATMAP_SCHEMA`. See [../export/formats.md](../export/formats.md) for the on-disk Parquet schemas |
| `export.py` | Versioned export dirs + shared orchestration; both frontends call `run_export()` so they produce byte-identical outputs. `ExportSelection` carries the request; `find_or_create_export_directory`/`commit_export_directory` do the SHA-256 settings-fingerprint versioning (reuse on match, else `_v2`, `_v3`, …); `RUN_SETTINGS_FIELDS` + `build_run_settings` serialise the run; `sanitize_filename` cleans names. See [../export/layout.md](../export/layout.md) |

`STATE_PREFIX` (the Streamlit session-state namespace) lives in
`streamlit_app/utils/state.py` — it's a UI-layer concern, not here.

---

## `core/spectra/`

mzML streaming + Koina AI prediction. Both produce the same `LoadedSpectrum`
so the annotation layer needn't know the spectrum's origin.

| Module | What it holds |
|---|---|
| `streaming_types.py` | `LoadedSpectrum` dataclass + the `StreamingResult` / `LoadResult` containers passed between stages |
| `spectrum_streaming.py` | `load_spectra_streaming()` — the public entry: builds an mzML → matches map, loads references first (so queries can be paired during streaming), then yields query results lazily |
| `mzml_native_ids.py` | Native-ID normalisation (`scan=`, `index=`, `spectrum=`). **timsTOF/Bruker `index=` IDs are 1-based** — the loader auto-detects the offset and shifts |
| `koina_prediction.py` | `predict_koina_spectrum()`, the `KOINA_MODELS` catalogue, `get_koina_settings()`, `convert_peptide_to_unimod()` |

---

## `core/plotting/`

One mirror renderer (Plotly + Matplotlib), two summary plot types, and the
parallel plot worker.

| Module | What it holds |
|---|---|
| `mirror.py` | Hand-rolled mirror renderer branching on `ms_matplotlib` vs `ms_plotly`. Pyopenms-viz' `plot_spectrum` was abandoned because it scrambles ion-type colouring when ion types are interleaved |
| `plot_worker.py` | `process_single_task()` — runs in a worker process: annotates the query/reference pair, renders, writes raster/vector/HTML. Carries the picklable `ComparisonTask` / `ComparisonResult` work units; format-specific save kwargs in `_FORMAT_SAVE_KWARGS` |
| `summary.py` | `generate_all_summary_plots()` over the `_PLOT_REGISTRY`, plus `plot_classic_dual_axis()` (PSM-counts bar + metric box-plot) and `plot_ridgeline()` (stacked per-peptidoform KDE). `METRIC_COLUMNS` selects which metrics summarise |

A "heatmap" summary plot is **not** generated — only `dual_axis` and
`ridgeline`. The openms_insight Heatmap on the Explore tab is a separate
interactive widget driven by `HEATMAP_SCHEMA`.

---

## `core/pipeline/`

The glue between lookup, load, annotate, compare, and plot. **This is the
layer the CLI and Streamlit share** — any feature that should run identically
from both front ends belongs here.

| Module | What it holds |
|---|---|
| `load_annotate.py` | The Load & Annotate driver. `run_streaming_annotation()` loops over `load_spectra_streaming` and applies `process_query_result()` (per streamed spectrum: annotate query, compute metrics, optional Koina fallback). `build_plot_tasks()` is the post-pass turning annotated queries into a `ComparisonTask` list |
| `kwargs_builders.py` | `build_annotation_kwargs()` / `build_koina_kwargs()` — take a generic `getter(key, default)` callable so the CLI (`SimpleNamespace`) and Streamlit (`params` dict) supply parameters without core knowing the storage format |
