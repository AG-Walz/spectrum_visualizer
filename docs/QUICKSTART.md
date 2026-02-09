# 🚀 SpecViz Quickstart Guide

Welcome to **SpecViz**, an MS/MS spectrum visualization tool for comparing experimental spectra with AI-powered predictions.

## Quick Start

### Installation

```bash
# Create conda environment
conda create -n specviz environment.yml
conda activate specviz
```

### Running the App

```bash
cd specviz
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

---

## 🧭 Navigate the App

The sidebar on the left contains all navigation options:

### 📄 Home Pages
- **Quickstart** (this page): Quick reference guide
- **Documentation**: Comprehensive guide with detailed explanations

### 🔬 Workspaces
Interactive environments for analysis:
- **Spectrum Analysis**: Main workspace for MS/MS visualization and comparison
- **Filter Large Datasets**: Preprocess mzML files before analysis
- **Settings**: Configure application preferences

**⚠️ Important**: Workspaces maintain state only during your session, unlike pages. Complete your analysis in one session or export results before navigating away, otherwise your progress will be overwritten.

---

## 🔬 Spectrum Analysis Workspace

The primary workspace for visualizing and analyzing MS/MS spectra.

### What it does:
- Loads idXML files containing peptide identifications
- Extracts unique peptidoforms and associated spectra
- Matches spectra to mzML files
- Generates AI predictions using Koina API (MS2PIP, Prosit, AlphaPept)
- Creates mirror plots comparing experimental vs predicted spectra
- Performs multi-spectrum comparisons (all-vs-all, one-vs-all)
- Exports plots and data to Excel/HTML

### When to use:
- Validating peptide identifications
- Comparing fragmentation patterns across replicates
- Evaluating AI model performance
- Quality control of MS/MS data

### Basic Workflow:

```
1. Upload idXML files
2. Select peptides of interest
3. Load associated spectra from mzML
4. Choose AI model for predictions
5. Generate comparisons (Exp vs Pred, All-vs-All)
6. Export results to Excel/HTML
```

---

## 🗂️ Filter Large Datasets Workspace

Preprocessing workspace for extracting relevant spectra from large mzML files.

### What it does:
- Filters mzML files based on peptide lists from idXML
- Extracts only relevant spectra (reduces file size 90-99%)
- Processes multiple files in batch
- Memory-efficient (handles 10GB files with 200MB RAM)

### When to use:
- Working with large mzML files (>1GB)
- Extracting spectra for specific peptides
- Reducing file sizes before sharing
- Preparing data for focused analysis

### Workflow:

```
1. Select idXML files
2. Define peptide filter
3. Process (auto-discovers mzML files)
4. Export filtered mzML files
```

---

## 📁 File Formats

### Identification Files (Required)
- **idXML** - OpenMS identification format (recommended)
- **mzIdentML** - PSI standard format
- Contains: Peptide sequences, modifications, spectrum references

### Spectrum Files (Required)
- **mzML** - Open mass spectrometry data format
- Contains: Raw MS/MS spectra (m/z, intensity pairs)

### Example File Structure:
```
data/
├── sample1_ids.idXML      # Peptide identifications
├── sample1_spectra.mzML   # Raw spectra
├── sample2_ids.idXML
└── sample2_spectra.mzML
```

---

## 🤖 AI Models

SpecViz uses the Koina API to access state-of-the-art fragmentation prediction models:

### MS2PIP (Recommended)
- **Best for**: timsTOF instruments, general use
- **Model**: ms2pip_timsTOF2024
- **Speed**: Fast (~1s per spectrum)
- **Accuracy**: High for modern instruments

### Prosit
- **Best for**: Orbitrap instruments
- **Model**: Prosit_2019_intensity
- **Speed**: Medium (~2s per spectrum)
- **Accuracy**: High for HCD fragmentation

### AlphaPept
- **Best for**: Generic predictions
- **Model**: AlphaPept_ms2_generic
- **Speed**: Medium (~2s per spectrum)
- **Accuracy**: Good across instruments

---

## 💡 Tips & Best Practices

### File Selection
- Use server browser for files >500MB
- Match idXML and mzML from same experiment
- Keep files in organized directory structure

### Peptide Selection
- Start with 5-10 peptides for testing
- Use "Select All" for comprehensive analysis
- Filter by modification if needed

### Analysis
- Choose AI model matching your instrument
- Use default settings for first analysis
- Adjust fragment tolerance if poor matches

### Export
- Export Excel for detailed data analysis
- Export HTML for sharing interactive plots
- Save results before changing workspaces

---

## 🆘 Troubleshooting

### "No matching mzML files found"
- Check file naming (should contain same base name)
- Verify files are in same directory
- Use server browser instead of upload

### "Koina API timeout"
- API may be busy, try again
- Reduce number of spectra
- Switch to different model

### "Memory error"
- Use Dataset Filter to reduce file sizes
- Process fewer files at once
- Close other applications

### "No spectra loaded"
- Check idXML contains spectrum references
- Verify mzML files are not corrupted
- Ensure matching run names

---

## 📚 Further Reading

For comprehensive documentation including:
- Advanced features and customization
- File format specifications
- API reference and technical details
- Performance optimization tips

See the **Documentation** page in the sidebar.

---

## 🔗 References

**SpecViz** is built on the OpenMS WebApps template:

Müller, T. D., Siraj, A., et al. (2025). OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. *Journal of Proteome Research*. [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)
