# MSInsight

MSInsight is a Streamlit application for interactive peptide-spectrum match (PSM)
analysis and visualization. It combines the openms_insight library for linked
interactive exploration with Plotly mirror plots, quantitative comparison metrics,
and batch export capabilities.

MSInsight replaces the earlier MSpecViz tool with a streamlined 4-tab workflow,
Koina AI spectrum prediction, and parallel export.

---

## Key Features

- **Interactive Explorer** -- openms_insight Table, LinePlot, SequenceView, and
  Heatmap components linked through a shared state manager. Click a PSM row to
  instantly view its annotated spectrum and sequence coverage.
- **Mirror Plots** -- Plotly interactive and Matplotlib static mirror plots
  comparing query spectra against experimental or AI-predicted references.
- **7 Comparison Metrics** -- R-squared (matching and binned), cosine similarity,
  spectral entropy similarity, Dice coefficient, Jaccard index, and shared peak
  ratio.
- **17 Ion Match Statistics** -- PyOpenMS SpectrumAnnotator metrics including
  coverage, annotated ion fraction, fragment mass errors, and series analysis.
- **Koina AI Prediction** -- 16 deep-learning models (MS2PIP, Prosit, AlphaPept,
  UniSpec, PredFull) for generating predicted reference spectra when no
  experimental reference is available.
- **Parallel Export** -- ProcessPoolExecutor-based batch generation of PNG, SVG,
  PDF, and HTML plots with versioned export directories and ZIP download.
- **Summary Plots** -- Dual-axis, ridgeline, and heatmap visualizations of
  metric distributions across peptidoforms.

---

## Prerequisites

- Python 3.11
- Conda (Miniconda or Anaconda)
- The `msinsight` conda environment

---

## Installation and Setup

Clone or copy the project, then activate the conda environment:

```bash
conda activate msinsight
```

The environment is located at `/mnt/volume/bin/miniconda/envs/msinsight/`.

Key dependencies (installed in the environment):

- streamlit
- openms_insight 0.1.11
- pyopenms
- polars
- plotly
- matplotlib
- ms_entropy
- koinapy
- psm_utils
- scipy

---

## Quick Start

Run the application:

```bash
streamlit run app.py
```

Or with the full path to the environment's Streamlit:

```bash
/mnt/volume/bin/miniconda/envs/msinsight/bin/streamlit run app.py
```

Open http://localhost:8501 in your browser. If connecting to a remote server,
set up an SSH tunnel:

```bash
ssh -L 8501:localhost:8501 user@server
```

See [docs/QUICKSTART.md](docs/QUICKSTART.md) for a walkthrough of the 4-tab
workflow with test data.

---

## Running Tests

```bash
# End-to-end pipeline test
python test_pipeline.py

# GUI smoke test
python test_gui.py

# Unit tests (if pytest is available)
pytest tests/
```

---

## Project Structure

| File / Directory | Description |
|---|---|
| `app.py` | Streamlit entry point |
| `content/peptide_mirror.py` | Main page: 4-tab workflow (Setup, Search & Load, Explore, Export) |
| `src/parquet_builder.py` | Converts PSMs and spectra to 4 Polars LazyFrames for openms_insight |
| `src/peptide_search.py` | Searches idXML files for target peptides |
| `src/spectrum_loader.py` | Streams spectra from mzML files |
| `src/peak_annotation.py` | Annotates spectra with theoretical fragment ions |
| `src/mirror_plot.py` | Plotly interactive and Matplotlib static mirror/single plots |
| `src/spectrum_metrics.py` | 7 comparison metrics and formatting |
| `src/summary_plots.py` | Dual-axis, ridgeline, and heatmap summary plots |
| `src/koina_prediction.py` | Koina AI spectrum prediction (16 models) |
| `src/parallel_processor.py` | Parallel plot generation and 17 ion match statistics |
| `src/export.py` | Versioned export directories with fingerprinting |
| `src/path_utils.py` | Path utility functions |
| `src/common/` | Shared Streamlit utilities |
| `utils/helpers.py` | Peptide list parsing, time formatting, nested dict access |
| `settings.json` | Application configuration (ion types, tolerances, Koina, export) |
| `test_pipeline.py` | End-to-end pipeline test |
| `test_gui.py` | GUI smoke test |

---

## Configuration

Application defaults are defined in `settings.json`:

- **ion_types** -- Default fragment ion types for annotation (default: b, y)
- **mz_tolerance** -- m/z matching tolerance in Da (default: 0.02)
- **annotation** -- Neutral losses, first prefix ion, precursor peak removal
- **koina** -- Default AI prediction model, collision energy, instrument type
- **performance** -- Worker count, cache size, chunk size
- **export** -- Default output directory, image formats, DPI, figure size
- **ion_colors** -- Color assignments for each ion type in plots

All settings can be overridden through the Setup tab in the UI.

---

## Documentation

| Guide | Description |
|---|---|
| [Quick Start](docs/QUICKSTART.md) | 4-tab workflow walkthrough with test data |
| [Metrics Reference](docs/METRICS.md) | All 7 comparison metrics and 17 IMS statistics with formulas |
| [Export Reference](docs/EXPORT.md) | Export directory structure, file formats, and schemas |
