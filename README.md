# MSpecViz

MSpecViz is a Streamlit-based web application for mass spectrometry (MS/MS) spectrum visualization and analysis. Compare experimental spectra against reference spectra with theoretical ion annotation.

Built on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/) [1].

## Table of Contents

- [Installation](#installation)
- [Basic Workflow](#basic-workflow)
- [Fast Loading with Parquet Conversion](#fast-loading-with-parquet-conversion)
- [Project Structure](#project-structure)
- [Testing](#testing)
- [Test Data](#test-data)
- [Documentation](#documentation)
  - [Caching System](docs/caching.md) — How auto-caching works, interrupted runs, invalidation
  - [Spectrum Comparison Metrics](docs/spectrum_metrics.md) — Cosine similarity, spectral entropy, R2, shared peaks
  - [Full User Guide](docs/DOCUMENTATION.md) — Comprehensive feature reference
  - [Quick Start](docs/QUICKSTART.md) — Get started in 5 minutes
  - [Parquet Conversion Tool](tools/README.md) — CLI reference for mzML-to-Parquet conversion
- [References](#references)

## Installation

### Requirements
- Python 3.11+
- Conda (Miniconda or Anaconda)

### Setup
```bash
# Create conda environment (includes all dependencies)
conda env create -f environment.yml
conda activate mspecviz

# Run the app
streamlit run app.py
```

Access at http://localhost:8501

### Remote Access (SSH tunnel)

If the server is behind a firewall (e.g., accessible only via SSH through a jumphost), direct access to port 8501 will be refused. Use an SSH tunnel to forward the port to your local machine.

**One command from your local PC:**

```bash
ssh -L 8501:localhost:8501 -J jumphost peptimine_highmem1
```

Then open http://localhost:8501 in your browser.

**Convenience script** — run `connect_mspecviz.sh` from your local machine to open the tunnel and launch the browser automatically:

```bash
# Copy scripts/connect_mspecviz.sh to your local machine, then:
chmod +x connect_mspecviz.sh
./connect_mspecviz.sh
```

> **Tip:** If port 8501 is already in use locally, pick another port:
> `ssh -L 9501:localhost:8501 -J jumphost peptimine_highmem1` then open http://localhost:9501

## Basic Workflow

1. **Files** — Select idXML (identifications) and mzML (spectra) files
2. **Peptides** — Enter target peptide sequences
3. **Reference** — Identify reference files by substring (e.g., "Synpep")
4. **Settings** — Configure ion types, m/z tolerance
5. **Search** — Find matching PSMs in idXML files
6. **Process** — Load spectra, annotate peaks, generate mirror plots, export

## Fast Loading with Parquet Conversion

For large datasets, convert your data to Parquet format **once** and load it instantly on every subsequent run. This skips all mzML parsing and index building.

### Option A: Command-line conversion (recommended for large datasets)

```bash
# Convert matched spectra to Parquet
conda run -n mspecviz python tools/convert_mzml_to_parquet.py \
  --idxml /path/to/*.idXML \
  --mzml-dir /path/to/mzml_files/ \
  --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
  --output exports/my_spectra.parquet \
  --reference-substring "Synpep"
```

Then in the app: go to **Files** tab, paste the `.parquet` path in the "Pre-converted Parquet" section, and jump directly to **Process**.

See [tools/README.md](tools/README.md) for full CLI documentation.

### Option B: Convert from within the app

1. Complete Steps 1-5 (Files through Search) as normal
2. In the Search results, click **"Convert to Parquet"**
3. On future runs, select the `.parquet` file in Step 1 and skip to Step 6

### Performance comparison

| Scenario | mzML (first run) | mzML (optimized index) | Parquet |
|----------|-------------------|------------------------|---------|
| 10 spectra, 1 file (100k spectra) | 5-30s | <1s | ~0.1s |
| 100 spectra, 7 files | 35-210s | ~4s | ~0.5s |
| Re-run same analysis | Same as first | Same as first | ~0.5s |

**Note:** After changing any code, restart the Streamlit app (`Ctrl+C` then `streamlit run app.py`) for changes to take effect. Streamlit caches modules in memory.

## Project Structure

```
MSpecViz/
├── app.py                          # Main entry point (st.navigation)
├── settings.json                   # App configuration
├── environment.yml                 # Conda dependencies
├── content/                        # Streamlit pages
│   ├── peptide_mirror.py               # Main analysis workflow
│   └── about.py                        # About page
├── src/                            # Backend processing
│   ├── common/common.py                # page_setup(), sidebar, UI helpers
│   ├── peptide_search.py               # idXML search
│   ├── spectrum_loader.py              # mzML + Parquet loading
│   ├── peak_annotation.py              # Ion annotation
│   ├── mirror_plot.py                  # Mirror plot generation
│   └── export.py                       # PNG/SVG/Excel export
├── tools/                          # Standalone scripts
│   └── convert_mzml_to_parquet.py      # mzML→Parquet converter
├── tests/                          # Unit tests
│   └── test_parquet_loader.py          # Parquet pipeline tests
├── utils/                          # Shared utilities
│   └── shared_imports.py               # Centralized heavy imports
└── docs/                           # Documentation
```

## Testing

```bash
# All tests
conda run -n mspecviz python -m pytest tests/ test_gui.py -v

# Parquet pipeline tests only
conda run -n mspecviz python -m pytest tests/test_parquet_loader.py -v

# GUI page tests only
conda run -n mspecviz python -m pytest test_gui.py -v
```

## Test Data

Location: `/mnt/volume/nacho/data/ProCa_Anna/small_dataset/mini_dataset`
- 7 mzML files (~2GB each)
- 3 idXML files
- Reference files contain "Synpep" in filename

## Documentation

| Guide | Description |
|-------|-------------|
| [Caching System](docs/caching.md) | Auto-caching, fingerprinting, interrupted runs, cache invalidation |
| [Spectrum Metrics](docs/spectrum_metrics.md) | Cosine similarity, spectral entropy, R2, shared peak ratios |
| [Full User Guide](docs/DOCUMENTATION.md) | Comprehensive reference for every feature and option |
| [Quick Start](docs/QUICKSTART.md) | Get up and running in 5 minutes |
| [Parquet Conversion](tools/README.md) | CLI tool for mzML-to-Parquet batch conversion |

## References

[1] Muller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. Journal of Proteome Research (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)
