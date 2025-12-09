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
├── .github/                     # GitHub workflows and templates
│   ├── ISSUE_TEMPLATE/          # Issue templates
│   └── workflows/               # CI/CD workflows
├── .streamlit/
│   └── config.toml              # Streamlit configuration
├── app.py                       # Main application entry point
├── Welcome.py                   # Welcome page
├── assets/                      # Static assets
├── cli/                         # Command-line tools
│   ├── filter_mzml.py           # mzML filtering CLI
│   ├── run_mzml_filter.sh       # Filter launch script
│   └── validate_filter.py       # Validation script
├── core/                        # Core processing modules
│   ├── correlation.py           # Correlation calculations
│   ├── mzml_processor.py        # mzML file processing
│   ├── spectrum_engine.py       # Spectrum analysis engine
│   └── streamlit_bridge.py      # Streamlit integration
├── data/                        # Data files and placeholders
├── docs/                        # Documentation
│   ├── QUICKSTART.md            # Quick start guide
│   └── DOCUMENTATION.md         # Full documentation
├── examples/                    # Example files
│   └── peptides.txt             # Example peptide list
├── pages/                       # Streamlit pages
│   ├── 1_Spectrum_Viewer.py     # Spectrum viewer page
│   ├── 2_Batch_Viewer.py        # Batch viewer page
│   ├── documentation.py         # Documentation page
│   └── quickstart.py            # Quickstart page
├── plots/                       # Plot outputs
├── utils/                       # Utility modules
│   ├── file_selection_components.py  # File selection UI
│   ├── helpers.py               # Helper functions
│   ├── shared_imports.py        # Shared imports
│   └── workspace_utils.py       # Workspace utilities
├── workspaces/                  # Analysis workspaces
│   ├── dataset_filter.py        # Dataset filtering
│   ├── settings.py              # Settings workspace
│   └── spectrum_analysis.py     # Spectrum analysis
├── environment.yml              # Conda environment
├── docker-compose.yml           # Docker configuration
└── README.md                    # This file
```


## References

[1] Müller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. Journal of Proteome Research (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)