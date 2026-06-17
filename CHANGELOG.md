# MSpecViz: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2026-06-17]

Initial release of MSpecViz is a peptide-spectrum match (PSM) analysis and visualization toolkit.
A pure-Python `core/` backend is shared by two interchangeable front ends: an
interactive Streamlit UI and a YAML-driven command line interface (CLI).

### `Added`

#### Architecture

- Pure-Python `core/` package with no front-end dependencies, consumed by both
  `streamlit_app/` and `cli/`.
- psm_utils `PSM`/`PSMList`-native data model streamed end-to-end.
- spectrum_utils ProForma annotation engine.
- Polars `LazyFrame`s feeding `openms_insight` components.
- Single-source `__version__` in `core/__init__.py`.

#### Peptide lookup & spectrum loading

- idXML peptide lookup with per-peptide PSM counting.
- Streaming mzML loader with native-ID resolution and `SpectrumLookup` fallback
  (handles 1-based Bruker/timsTOF native IDs).
- Best-spectrum selection by precursor m/z and retention-time proximity.
- Recursive directory scanning for idXML and mzML inputs.
- Parallel loading/search across all pipeline stages via joblib.

#### Annotation & metrics

- Theoretical fragment generation with internal-fragment ion support.
- 7 quantitative comparison metrics (cosine, R², spectral entropy, …) plus
  ion-match statistics.

#### Koina AI prediction

- Koina spectrum prediction (14 models) as fallback reference for unpaired
  spectra, with configurable model, collision energy, and instrument.

#### Visualization

- Plotly interactive mirror / single-spectrum plots via pyopenms-viz.
- Matplotlib static plots for batch export.
- Summary plots: dual-axis, ridgeline (seaborn/joypy), heatmap.
- Per-peak coloring from a shared `ION_COLORS` palette.

#### Export

- Versioned export directories with `run_settings.json` fingerprinting.
- Output formats: Parquet, CSV/XLSX, PNG/SVG/PDF, and interactive HTML.
- Parallel plot generation with RAM-backed tmpfs staging.
- Resumable Parquet caching at multiple pipeline stages.
- Plot gallery browser (grouped by peptidoform → reference), export-history
  browser, and ZIP download.

#### Streamlit UI (4-tab workflow)

- **Setup**: directory file browser, peptide picker with PSM count and min-PSM
  filter, annotation/Koina/performance settings.
- **Load & Annotate**: lookup, streaming annotation, metric computation.
- **Explore**: linked Interactive Explorer (Table + SequenceView + LinePlot +
  Heatmap), peptide-scoped Mirror Plots picker, Comparison Metrics view.
- **Export**: versioned exports, gallery, and history.
- Workspace persistence via `params.json`; `@st.fragment` isolation to minimize
  reruns; multi-instance support.

#### CLI

- Headless pipeline driven by a single YAML config (`config_example.yaml`).
- click-based interface with tqdm progress bars.
- idXML/mzML filter flags, `--no-plots`, and `--workers`.

#### Filtering

- `--idxml-filter` / `--mzml-filter` (and UI equivalents) matching against full
  paths to subset inputs.

#### Testing & CI

- pytest suite nested under `tests/core`, `tests/streamlit`, `tests/cli`,
  including GUI (`AppTest`) and pipeline characterization tests.
- GitHub Actions CI workflow.

#### Documentation

- `docs/` set (architecture, cli, export, metrics, settings) and a README with
  table of contents.

### `Dependencies`

| Dependency       | Note                                          |
| ---------------- | --------------------------------------------- |
| Python           | 3.11                                          |
| openms_insight   | 0.1.13 — linked Table/SequenceView/LinePlot   |
| streamlit        | >=1.48.0 — front-end UI                        |
| polars           | LazyFrame data model                          |
| psm_utils        | native PSM/PSMList model                       |
| spectrum_utils   | ProForma annotation                           |
| pyopenms         | spectrum I/O and lookup                        |
| pyopenms-viz     | Plotly mirror-plot rendering                   |
| koina            | AI spectrum prediction                        |
| joblib           | parallel pipeline stages                       |
| click / tqdm     | CLI interface and progress                    |
| seaborn / joypy  | summary plots                                 |
| ms_entropy       | spectral-entropy metric                       |
