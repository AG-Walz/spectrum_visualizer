# SpecViz App

SpecViz App is designed to help you analyze and visualize your mass spectrometry data. It provides a user-friendly interface for data upload, filtering, and visualization.

This app is built based on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/) [1].

## Installation

### Requirements
Ensure python \leq 3.11+ is installed. Sufficient RAM to load the largest file. Required packages are listed under [environment.yml](environment.yml)

### Setup
```bash
# Create conda environment
conda env create -f environment.yml
conda activate specviz

# Install dependencies
pip install pyopenms streamlit psm_utils pandas numpy matplotlib openpyxl
```

## Quick Start

### Spectrum Visualization
```bash
cd specviz
streamlit run app.py
```
Access at http://localhost:8501

## Documentation

- **Quickstart Guide**: [docs/QUICKSTART.md](docs/QUICKSTART.md) - Get started quickly
- **Full Documentation**: [docs/DOCUMENTATION.md](docs/DOCUMENTATION.md) - Comprehensive guide

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