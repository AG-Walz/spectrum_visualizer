# SpecViz App

SpecViz App is designed to help you analyze and visualize your mass spectrometry data. It provides a user-friendly interface for data upload, filtering, and visualization.

This app is built based on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/) [1].

## Features

### 1. Spectrum Visualization Workspace
- **Spectral alignment analysis**: Compare MS/MS spectra across samples
- **Peptide comparison modes**: One-vs-One, All-vs-All, One-vs-All
- **R-squared metrics**: Matching ions and binned intensity correlation
- **Export functionality**: Excel reports with overview plots
- **Smart file matching**: Automatic PSM to mzML file association
- **Progress tracking**: Real-time feedback for uploads and validation

### 2. mzML Filter Workspace (NEW)
- **Peptide-based filtering**: Create compact mzML files with only relevant spectra
- **Batch processing**: Process multiple idXML and mzML files automatically
- **Modified peptide matching**: Automatically includes all modifications
- **File consolidation**: Append spectra across runs into single files
- **Web interface and CLI**: Choose your preferred workflow

## Quick Start

### Spectrum Visualization
```bash
cd specviz
streamlit run app.py
```
Access at http://localhost:8501

### mzML Filtering

**Web Interface:**
```bash
cd specviz
./cli/run_mzml_filter.sh
```
Access at http://localhost:8502

**Command Line:**
```bash
cd specviz
python cli/filter_mzml.py \
    -i path/to/idxml \
    -m path/to/mzml \
    -o path/to/output \
    -p PEPTIDEK,SEQUENCEM,PROTEINR
```

## Documentation

- **Quickstart Guide**: [docs/QUICKSTART.md](docs/QUICKSTART.md) - Get started quickly
- **Full Documentation**: [docs/DOCUMENTATION.md](docs/DOCUMENTATION.md) - Comprehensive guide

## Installation

### Requirements
- Python 3.8+
- PyOpenMS 3.4.0+
- Streamlit
- psm_utils
- pandas, numpy, matplotlib

### Setup
```bash
# Create conda environment
conda create -n specviz python=3.11
conda activate specviz

# Install dependencies
pip install pyopenms streamlit psm_utils pandas numpy matplotlib openpyxl
```

## Project Structure

```
specviz/
├── bin/
│   ├── SpecViz.py              # Main spectrum visualization app
│   ├── mzml_filter.py          # mzML filtering module
│   ├── mzml_filter_cli.py      # CLI for mzML filtering
│   └── run_mzml_filter.sh      # Launch script for filter workspace
├── src/
│   ├── spectrum_workspace.py   # Spectrum analysis workspace
│   └── mzml_filter_workspace.py # mzML filtering workspace
├── docs/
│   └── MZML_FILTER.md          # mzML filter documentation
└── README.md                    # This file
```


## References

[1] Müller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. Journal of Proteome Research (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)