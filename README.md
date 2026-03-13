# MSpecViz

MSpecViz is a Streamlit-based web application for mass spectrometry (MS/MS) spectrum visualization and analysis. Compare experimental spectra against reference spectra with theoretical ion annotation, generate publication-ready mirror plots, and compute quantitative similarity metrics.

Built on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/) [1].

---

## Table of Contents

- **Getting Started**
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Two ways to run MSpecViz](#two-ways-to-run-mspecviz) — Streamlit app and CLI pipeline
  - [Remote Access (SSH Tunnel)](#remote-access-ssh-tunnel)
  - [Quick Start Guide](docs/QUICKSTART.md) -- 5-minute walkthrough
- **Usage Modes**
  - [Streamlit Web App](#streamlit-web-app) -- Interactive browser-based analysis
  - [CLI Pipeline](#cli-pipeline) -- Terminal-based batch processing
  - [Docker](#docker) -- Containerized deployment
- **Workflow & Configuration**
  - [Input Files](#input-files) -- What you need and in what format
  - [Analysis Decisions](#analysis-decisions) -- Choices you need to make before running
  - [Fast Loading with Parquet](#fast-loading-with-parquet) -- Speed up repeated analyses 100x
- **Detailed Guides**
  - [Streamlit Workflow Steps](#streamlit-workflow-steps)
  - [CLI Configuration](#cli-configuration)
  - [Comparison Metrics](#comparison-metrics)
- **Standalone Scripts**
  - [Summary & Per-Peptide Plots](#summary--per-peptide-plots) -- Cross-peptidoform summaries and per-peptide dumbbell plots
  - [Parquet Conversion Tool](#parquet-conversion-tool) -- Batch mzML-to-Parquet conversion
- **Reference**
  - [Full User Guide](docs/DOCUMENTATION.md) -- Comprehensive reference for every feature
  - [Detailed Documentation Index](#detailed-documentation-index)
  - [Testing](#testing)

---

## Requirements

- **Conda** (Miniconda or Anaconda) — manages all dependencies
- **Python 3.11** — pinned by `environment.yml`
- Core dependencies (installed automatically): PyOpenMS, Streamlit ≥ 1.48, Pandas, NumPy, Bokeh, Plotly, psm_utils, scikit-learn, koinapy

No manual package installation is needed — `environment.yml` pins everything.

## Installation

```bash
git clone <repository-url>
cd MSpecViz
conda env create -f environment.yml   # creates the 'mspecviz' environment
```

This downloads and installs all dependencies into an isolated `mspecviz` conda environment. Subsequent commands must be run inside that environment (see below).

## Two ways to run MSpecViz

MSpecViz has two interfaces that share the same analysis backend — results are identical regardless of which you use:

| Interface | Best for |
|-----------|----------|
| **Streamlit web app** | Interactive analysis, browsing results in real time |
| **CLI pipeline** | Batch processing, scripting, HPC job submission |

### Streamlit Web App

```bash
conda activate mspecviz
streamlit run app.py
```

Open http://localhost:8501 in your browser. For a full walkthrough see the [Quick Start Guide](docs/QUICKSTART.md).

Alternatively, without activating the environment first:
```bash
conda run -n mspecviz streamlit run app.py
```

### CLI Pipeline

```bash
conda activate mspecviz
python cli.py --config examples/config_cli.yaml
```

The CLI reads a YAML config file and runs the full pipeline without a browser. See [CLI Pipeline](#cli-pipeline) below or the [CLI Reference](docs/cli.md) for the complete documentation.

Alternatively, without activating first:
```bash
conda run -n mspecviz python cli.py --config examples/config_cli.yaml
```

### Remote Access (SSH Tunnel)

If the server is behind a firewall (accessible only via SSH through a jumphost):

```bash
ssh -L 8501:localhost:8501 -J jumphost server_name
```

Then open http://localhost:8501 in your local browser.

> **Tip:** If port 8501 is already in use locally, pick another port:
> `ssh -L 9501:localhost:8501 -J jumphost server_name` then open http://localhost:9501

A convenience script is provided at `scripts/connect_mspecviz.sh`.

---

## Input Files

MSpecViz requires two types of input files from the same experiment:

| File Type | Extension | Contains | Created By |
|-----------|-----------|----------|------------|
| **Identification files** | `.idXML` | Peptide-spectrum matches (PSMs) from a search engine | OpenMS (MSGF+, Comet, etc.) |
| **Spectrum files** | `.mzML` | Raw MS/MS spectra (m/z + intensity arrays) | Instrument conversion (e.g., msconvert) |

The idXML files contain metadata linking each PSM to a spectrum inside a specific mzML file. Both must come from the same experiment.

**Optional**: A pre-converted `.parquet` file from a previous run or the [Parquet conversion tool](tools/README.md) can replace both inputs for instant loading.

---

## Analysis Decisions

Before running an analysis (either Streamlit or CLI), you need to decide on the following:

### 1. Target Peptides

Which peptide sequences to analyze. You can provide:
- A list of base sequences (unmodified) -- modified variants are found automatically
- A file with one peptide per line (see `examples/peptides.txt`)
- `all` (CLI only) -- extract every peptide found in the idXML files

### 2. Reference Identification

A **substring** that identifies your reference mzML files (e.g., `Synpep`, `Standard`, `Ref`). Files whose name contains this substring become **reference spectra**; all others become **queries**. Mirror plots pair each query with its matching reference.

If omitted, every spectrum is plotted individually with annotations (no mirror comparison). See [Reference Setup](docs/reference_setup.md).

### 3. Ion Types and Tolerance

| Setting | Default | Description | Guide |
|---------|---------|-------------|-------|
| **Ion types** | b, y | Which fragment ion series to annotate (a, b, c, x, y, z, precursor, immonium) | [Annotation Settings](docs/annotation_settings.md) |
| **m/z tolerance** | 0.02 Da | Maximum distance between experimental and theoretical peaks for a match | [Annotation Settings](docs/annotation_settings.md) |
| **Neutral losses** | Off | Include H2O and NH3 loss peaks | [Annotation Settings](docs/annotation_settings.md) |
| **Remove precursor** | Off | Strip precursor ions from plots | [Annotation Settings](docs/annotation_settings.md) |
| **Intensity mode** | Relative | Normalize to base peak (Relative) or keep raw values (Raw) | [Annotation Settings](docs/annotation_settings.md) |

### 4. Export Settings

| Setting | Default | Description |
|---------|---------|-------------|
| **Project name** | `MSpecViz_analysis` | Base name for the export folder |
| **Export directory** | `<workspace>/exports` | Where output is saved (defaults to the active workspace; see [Workspaces](docs/workspaces.md)) |
| **Save plots to disk** | On | Generate and save individual plot files; disable for summary-CSV-only runs (`--no-plots` in CLI) |
| **Formats** | PNG | Image formats: PNG (raster, 300 DPI) and/or SVG (vector); only visible when Save plots is on |
| **Plot unpaired** | Off | Whether to plot spectra that have no matching reference |

See [Export & Output](docs/export.md) for the full output directory structure.

---

## Streamlit Web App

The Streamlit interface provides an interactive 7-step workflow with real-time feedback, progress bars, and in-browser result preview.

### Streamlit Workflow Steps

| Step | Tab | What You Do | Details |
|------|-----|-------------|---------|
| 1 | **Files** | Select idXML and mzML files using the directory browser | [File Selection](docs/file_selection.md) |
| 2 | **Peptides** | Enter target peptide sequences (one per line) | [Peptide Search](docs/peptide_search.md) |
| 3 | **Reference** | Enter a substring to identify reference files | [Reference Setup](docs/reference_setup.md) |
| 4 | **Settings** | Configure ion types, m/z tolerance, losses | [Annotation Settings](docs/annotation_settings.md) |
| 5 | **Search** | Click "Search Peptides" to find matching PSMs | [Peptide Search](docs/peptide_search.md) |
| 6 | **Process** | Configure export and click "Process & Export" | [Export & Output](docs/export.md) |
| 7 | **Results** | View summary of exported files | [Mirror Plots](docs/mirror_plots.md) |

Each step validates its prerequisites before allowing you to proceed. You can go back to earlier tabs to adjust settings at any time.

**Shortcut**: If you have a pre-converted Parquet file, paste its path in Step 1 and skip directly to Step 6.

---

## CLI Pipeline

The CLI runs the same analysis from the terminal without opening a browser. Results are identical to the Streamlit app.

```bash
# Basic run with a config file
python cli.py --config examples/config_cli.yaml

# Verbose output for debugging
python cli.py --config examples/config_cli.yaml --verbose

# Override settings on the command line
python cli.py --config examples/config_cli.yaml --mz-tolerance 0.05 --formats png svg

# Speed up large analyses with more parallel workers
python cli.py --config examples/config_cli.yaml --workers 8
```

### CLI Configuration

Configuration is merged from three layers (highest priority wins):

| Priority | Source | Example |
|----------|--------|---------|
| 1 (highest) | CLI arguments | `--mz-tolerance 0.05` |
| 2 | YAML config file | `annotation: mz_tolerance: 0.02` |
| 3 (lowest) | `settings.json` | App-wide defaults |

Copy `examples/config_cli.yaml` to get started -- it contains every option with inline documentation.

### Common CLI Overrides

```bash
# Use a peptide list file
python cli.py --config my_config.yaml --peptides my_peptides.txt

# Process all peptides found in the idXML files
python cli.py --config my_config.yaml --peptides all

# Export both PNG and SVG
python cli.py --config my_config.yaml --formats png svg

# Speed up plot generation with more parallel workers (default: 4)
python cli.py --config my_config.yaml --workers 8

# Re-run from cached Parquet (instant)
python cli.py --config my_config.yaml --parquet exports/previous_run/.mspecviz_cache/abc123.parquet

# Narrow a large glob to a subset for debugging (case-insensitive filename match)
python cli.py --config my_config.yaml --idxml-filter run1 --mzml-filter run1

# Quiet mode for scripting (exit codes: 0=success, 1=failure, 2=config error)
python cli.py --config my_config.yaml --quiet
```

See the [CLI Reference](docs/cli.md) for the full list of arguments, pipeline phases, and output structure.

---

## Standalone Scripts

In addition to the main CLI pipeline (`cli.py`), MSpecViz includes two standalone scripts for post-processing and data conversion.

### Summary & Per-Peptide Plots

`src/summary_plots.py` generates visualizations from the `comparison_summary.csv` produced by the main pipeline.

**Summary plots** (cross-peptidoform): dual-axis bar+box, ridgeline distribution, and heatmap plots showing a chosen metric across all peptidoforms.

**Per-peptide dumbbell plots**: For each peptidoform, shows R² (binned) vs R² (match) as connected dots for every query, faceted by reference. Highlights where the two R² metrics diverge.

```bash
# Summary plots (all types)
python -m src.summary_plots exports/my_run/comparison_summary.csv

# Per-peptide dumbbell plots
python -m src.summary_plots exports/my_run/comparison_summary.csv --peptide-detail

# List available metrics
python -m src.summary_plots exports/my_run/comparison_summary.csv --list-metrics
```

Summary plots are saved to `{export_dir}/summary/`. Dumbbell plots are saved alongside mirror plots in `{export_dir}/png/{Peptide}/`.

### Parquet Conversion Tool

`tools/convert_mzml_to_parquet.py` converts idXML + mzML files into a single Parquet file for instant loading.

```bash
python tools/convert_mzml_to_parquet.py \
  --idxml /path/to/*.idXML \
  --mzml-dir /path/to/mzml_files/ \
  --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
  --output exports/my_spectra.parquet \
  --reference-substring "Synpep"
```

See [CLI Reference: Standalone Scripts](docs/cli.md#standalone-scripts) for the full list of arguments, and [Parquet Conversion Tool](tools/README.md) for detailed documentation.

---

## Docker

```bash
# Build and run
docker compose up --build

# Or just build
docker build -t mspecviz .
docker run -p 8501:8501 -v ./data:/app/data mspecviz
```

The container uses the conda environment defined in `environment.yml`.

---

## Fast Loading with Parquet

For large datasets, convert your data to Parquet format **once** and load it instantly on every subsequent run. This skips all mzML parsing and index building.

| Scenario | mzML (first run) | mzML (cached index) | Parquet |
|----------|-------------------|---------------------|---------|
| 10 spectra, 1 file (100k spectra) | 5-30s | <1s | ~0.1s |
| 100 spectra, 7 files | 35-210s | ~4s | ~0.5s |
| Re-run same analysis | Same as first | Same as first | ~0.5s |

### Three ways to use Parquet

1. **Auto-cache** (default): After a successful run, MSpecViz automatically saves a Parquet cache. The next run with the same inputs loads from cache instead of mzML. See [Caching System](docs/caching.md).

2. **In-app conversion**: After Step 5 (Search) in the Streamlit app, click "Convert to Parquet" to save matched spectra. Paste the Parquet path in Step 1 on future runs.

3. **CLI conversion tool**: For large datasets or batch processing:
   ```bash
   python tools/convert_mzml_to_parquet.py \
     --idxml /path/to/*.idXML \
     --mzml-dir /path/to/mzml_files/ \
     --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
     --output exports/my_spectra.parquet \
     --reference-substring "Synpep"
   ```
   See [Parquet Conversion Tool](tools/README.md) for full documentation.

---

## Comparison Metrics

Seven quantitative metrics are computed for each query-reference pair:

| Metric | Scope | What It Answers |
|--------|-------|-----------------|
| R2 Matching | Annotated ions | Do shared fragment ions have similar relative intensities? |
| Cosine Similarity | Annotated ions | How similar are the intensity vectors of shared ions? |
| Dice Coefficient | Annotated ions | What fraction of annotated ions are shared? |
| Jaccard Index | Annotated ions | Same as Dice but with stricter union denominator |
| Shared Peak Ratio | Annotated ions | Is the smaller spectrum's information contained in the larger? |
| R2 Binned | All peaks | Do the overall spectral profiles correlate? |
| Spectral Entropy | All peaks | Information-theoretic similarity across all peaks |

Metrics appear on each mirror plot and in the `comparison_summary.csv` / `.xlsx` export. See [Spectrum Metrics](docs/spectrum_metrics.md) for equations and interpretation.

---

## Testing

```bash
# All tests
conda run -n mspecviz python -m pytest tests/ test_gui.py -v

# Unit tests only
conda run -n mspecviz python -m pytest tests/ -v

# GUI page tests only
conda run -n mspecviz python -m pytest test_gui.py -v
```

---

## Detailed Documentation Index

### Getting Started

| Guide | Description |
|-------|-------------|
| [Quick Start](docs/QUICKSTART.md) | 5-minute walkthrough of the full workflow |
| [Full User Guide](docs/DOCUMENTATION.md) | Comprehensive reference for every feature, setting, and workflow step |
| [Workspaces](docs/workspaces.md) | What workspaces are, how to create/switch/delete them, sharing, and CLI integration |

### Workflow Steps (Streamlit)

| Guide | Description |
|-------|-------------|
| [File Selection](docs/file_selection.md) | Directory browsing, file disambiguation, Parquet shortcut, mzML resolution |
| [Peptide Search](docs/peptide_search.md) | How idXML search works, PSM extraction, variant grouping |
| [Reference Setup](docs/reference_setup.md) | Reference vs query categorization, unpaired spectra, choosing a substring |
| [Annotation Settings](docs/annotation_settings.md) | Ion types, m/z tolerance, neutral losses, precursor handling, intensity modes |
| [Mirror Plots](docs/mirror_plots.md) | Plot anatomy, how to interpret colors and labels, single spectrum plots |
| [Export & Output](docs/export.md) | Directory structure, plot formats, comparison summary spreadsheet |

### Technical References

| Guide | Description |
|-------|-------------|
| [CLI Reference](docs/cli.md) | Full CLI documentation: config options, arguments, pipeline phases, output structure, standalone scripts |
| [Spectrum Metrics](docs/spectrum_metrics.md) | Seven comparison metrics with equations: R2, cosine, entropy, Dice, Jaccard, shared ratio |
| [Caching System](docs/caching.md) | Auto-caching, fingerprinting, interrupted runs, cache invalidation |
| [Parquet Conversion Tool](tools/README.md) | CLI tool for standalone mzML-to-Parquet batch conversion |
| [Example CLI Config](examples/config_cli.yaml) | Annotated YAML configuration file with every available option |

### Standalone Scripts

| Script | Description |
|--------|-------------|
| [`cli.py`](docs/cli.md) | Main analysis pipeline: search, load, annotate, plot, export |
| [`src/summary_plots.py`](docs/cli.md#summary--per-peptide-plots-srcsummary_plotspy) | Summary plots (dual-axis, ridgeline, heatmap) and per-peptide dumbbell plots |
| [`tools/convert_mzml_to_parquet.py`](tools/README.md) | Batch mzML+idXML → Parquet conversion for instant loading |

---

## References

[1] Muller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. Journal of Proteome Research (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)
