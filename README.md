# MSpecViz

MSpecViz is a Streamlit-based web application for mass spectrometry (MS/MS) spectrum visualization and analysis. Compare experimental spectra against reference spectra with theoretical ion annotation, generate publication-ready mirror plots, and compute quantitative similarity metrics.

Built on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/) [1].

---

## Table of Contents

- **Getting Started**
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Running the App](#running-the-app)
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
- **Reference**
  - [Full User Guide](docs/DOCUMENTATION.md) -- Comprehensive reference for every feature
  - [Detailed Documentation Index](#detailed-documentation-index)
  - [Testing](#testing)

---

## Requirements

- Python 3.11+
- Conda (Miniconda or Anaconda)

## Installation

```bash
git clone <repository-url>
cd MSpecViz
conda env create -f environment.yml
```

## Running the App

### Streamlit Web App

```bash
conda activate mspecviz
streamlit run app.py
```

Open http://localhost:8501 in your browser.

### CLI Pipeline

```bash
conda activate mspecviz
python cli.py --config examples/config_cli.yaml
```

See [CLI Pipeline](#cli-pipeline) below or the [CLI Reference](docs/cli.md) for the complete documentation.

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
| **Export directory** | `./exports` | Where output is saved |
| **Formats** | PNG | Image formats: PNG (raster, 300 DPI) and/or SVG (vector) |
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

# Re-run from cached Parquet (instant)
python cli.py --config my_config.yaml --parquet exports/previous_run/.mspecviz_cache/abc123.parquet

# Quiet mode for scripting (exit codes: 0=success, 1=failure, 2=config error)
python cli.py --config my_config.yaml --quiet
```

See the [CLI Reference](docs/cli.md) for the full list of arguments, pipeline phases, and output structure.

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
| [CLI Reference](docs/cli.md) | Full CLI documentation: config options, arguments, pipeline phases, output structure |
| [Spectrum Metrics](docs/spectrum_metrics.md) | Seven comparison metrics with equations: R2, cosine, entropy, Dice, Jaccard, shared ratio |
| [Caching System](docs/caching.md) | Auto-caching, fingerprinting, interrupted runs, cache invalidation |
| [Parquet Conversion Tool](tools/README.md) | CLI tool for standalone mzML-to-Parquet batch conversion |
| [Example CLI Config](examples/config_cli.yaml) | Annotated YAML configuration file with every available option |

---

## References

[1] Muller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. Journal of Proteome Research (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)
