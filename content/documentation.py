import streamlit as st

# Simple setup - no need for complex common.py
st.set_page_config(
    page_title="SpecViz - Documentation (under construction)",
    initial_sidebar_state="expanded"
)

st.title("📖 SpecViz Documentation")

# Table of Contents
st.sidebar.markdown("## 📋 Table of Contents")
sections = [
    "Overview",
    "File Formats",
    "AI Models",
    "Workflow",
    "Features",
    "Troubleshooting",
    "Technical Details",
    "API Reference"
]

for section in sections:
    st.sidebar.markdown(f"- [{section}](#{section.lower().replace(' ', '-')})")

# Main documentation content
st.markdown("---")  # Separator
st.markdown("""
## Overview

SpecViz is a modern web application for visualizing MS/MS spectra with AI-powered predictions. It provides an intuitive interface for comparing experimental mass spectrometry data with theoretical predictions from state-of-the-art machine learning models.

---

## File Formats

### Supported Input Formats

**mzML Files**
- Open standard format for mass spectrometry data
- Contains MS/MS spectra with m/z values and intensities
- Must contain MS2 spectra for visualization
- Supports various instrument types (Orbitrap, timsTOF, etc.)

**Identification Files**
- **idXML**: OpenMS identification format
- **mzIdentML**: PSI standard for peptide identifications
- **mzid**: Alternative extension for mzIdentML

### Data Requirements

- Files should be from the same experiment
- ID files must contain peptide sequences with spectrum references
- Spectra must have corresponding identifications for visualization

---

## AI Models

SpecViz integrates with the Koina API to provide access to various AI models:

### MS2PiP Models
- **ms2pip_timsTOF2024**: Latest model optimized for timsTOF instruments
- **ms2pip_HCD2021**: HCD fragmentation model from 2021

### Prosit Models
- **Prosit_2019_intensity**: Original Prosit model for intensity prediction
- **Prosit_2020_intensity**: Updated model with improved accuracy

### AlphaPept Models
- **AlphaPept_ms2_generic**: Generic fragmentation model

---

## Workflow

### 1. Data Upload
Navigate to **Workspaces** → **Spectrum Analysis** and upload your files:
- Select mzML file containing MS/MS spectra
- Select idXML/mzIdentML file with peptide identifications

### 2. Model Selection
Choose an AI model from the sidebar:
- Default: ms2pip_timsTOF2024 (recommended for most data)
- Select based on your instrument type and data characteristics

### 3. Spectrum Selection
- Browse the interactive table of identified spectra
- View key properties: peptide sequence, charge, RT, m/z
- Select a spectrum index for visualization

### 4. Visualization
- Compare experimental (top) vs predicted (bottom) spectra
- Analyze ion annotations and intensities
- Use interactive features for detailed exploration

---

## Features

### Visualization
- **Mirror plots**: Side-by-side comparison of experimental vs predicted spectra
- **Ion annotation**: Automatic b-ion and y-ion identification
- **Color coding**: Distinct colors for different ion types
- **Interactive legends**: Toggle visibility of ion types

### Peptide Modifications
- Full support for modified peptides
- Proper mass shift formatting (e.g., `M[+15.9949]` for oxidation)
- Automatic conversion between formats

### Error Handling
- Comprehensive file validation
- Multiple visualization fallback methods
- Detailed error messages and troubleshooting tips

### Performance
- Efficient file caching
- Optimized data processing
- Responsive UI with progress indicators

## Troubleshooting

### Common Issues

**"No spectra found"**
- Ensure mzML file contains MS2 spectra
- Check that MS level is 2 in the data

**"No identifications found"**
- Verify ID file contains peptide sequences
- Check file format compatibility

**"No matching spectra"**
- Ensure files are from the same experiment
- Check spectrum indexing consistency

**Visualization errors**
- Enable debug mode for detailed information
- Check browser compatibility
- Try different visualization methods

### Debug Mode
Enable debug mode in the sidebar to:
- View detailed processing information
- See intermediate data structures
- Identify specific error locations

## Technical Details

### Architecture
- **Frontend**: Streamlit web framework
- **Backend**: OpenMS for MS data processing
- **AI Integration**: Koina API for model access
- **Visualization**: Bokeh for interactive plots

### Dependencies
- OpenMS: Mass spectrometry data processing
- Pandas: Data manipulation and analysis
- Bokeh: Interactive visualization
- Streamlit: Web application framework

### File Processing
1. **Validation**: Check file formats and content
2. **Parsing**: Extract spectra and identifications
3. **Mapping**: Link spectra to peptide identifications
4. **Normalization**: Standardize peptide sequences

### Visualization Pipeline
1. **Data extraction**: Get experimental spectrum
2. **Prediction**: Query AI model via Koina
3. **Processing**: Normalize and annotate spectra
4. **Rendering**: Generate interactive mirror plot

## API Reference

### Core Functions

**File I/O**
- `load_experiment(path)`: Load mzML file
- `load_identifications(path)`: Load ID file
- `extract_spectrum_mapping()`: Create spectrum-peptide mapping

**Data Processing**
- `normalize_peptide_sequence()`: Standardize peptide format
- `convert_to_mass_shift_format()`: Convert modification format
- `validate_files()`: Check file integrity

**Visualization**
- `visualize_spectrum()`: Create mirror plot
- `display_plot_with_fallbacks()`: Render with multiple methods
- `show_visualization_info()`: Display plot metadata

### Configuration

Settings are stored in `settings.json`:
```json
{
    "app-name": "SpecViz",
    "version": "1.0.0",
    "models": {
        "default": "ms2pip_timsTOF2024",
        "available": [...]
    },
    "ui": {
        "theme": "light",
        "debug_mode": false
    }
}
```

### Extension Points
- Add new AI models via Koina integration
- Extend file format support
- Customize visualization themes
- Add new analysis features

---

*For additional help, please check the GitHub repository or contact the development team.*
""")