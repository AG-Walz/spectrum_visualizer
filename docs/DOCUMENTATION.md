# 📖 SpecViz Documentation

*Comprehensive guide to every feature, option, and workflow in the application*

---

## Table of Contents

1. [Overview](#overview)
2. [File Upload](#file-upload)
3. [idXML Processing](#idxml-processing)
4. [Peptide Selection](#peptide-selection)
5. [mzML Loading](#mzml-loading)
6. [Analysis Configuration](#analysis-configuration)
7. [Automated Analysis](#automated-analysis)
8. [Visualization](#visualization)
9. [Export Options](#export-options)
10. [File Formats](#file-formats)
11. [AI Models](#ai-models)
12. [Troubleshooting](#troubleshooting)
13. [Advanced Features](#advanced-features)
14. [Technical Details](#technical-details)

---

## Overview

SpecViz is a streamlined web application for MS/MS spectrum visualization with AI-powered fragmentation predictions. It compares experimental mass spectrometry data against predictions from state-of-the-art machine learning models (MS2PIP, Prosit, AlphaPept) via the Koina API.

### Core Capabilities
- Load peptide identifications from idXML/mzIdentML files
- Extract and filter peptidoforms (sequences with modifications)
- Match spectra to mzML files efficiently
- Generate AI-based fragmentation predictions
- Create interactive mirror plots with ion annotations
- Perform multi-spectrum comparisons (replicate analysis)
- Export results to Excel and HTML formats

### Typical Workflow
```
Upload Files → Select Peptides → Load Spectra → Configure Analysis → Generate Plots → Export Results
```

---

## File Upload

The first step in any analysis is uploading your data. SpecViz uses a two-phase approach to optimize memory usage.

### Upload Methods

**Option 1: Browser Upload**
- Click "🔍 Upload Identification Files"
- Select files from your computer
- Best for: Small to medium files (<500 MB)
- Limitation: Browser memory limits may apply

**Option 2: Server Browser**
- Click "📂 Browse Server Files"
- Navigate to data directory
- Select files directly from server filesystem
- Best for: Large files, batch processing, shared datasets

### Supported File Formats

**Identification Files (Required):**
- `.idxml` - OpenMS identification format (recommended)
- `.mzid` / `.mzidentml` - PSI standard format
- Contains: Peptide sequences, modifications, scores, spectrum references

**Spectrum Files (Required):**
- `.mzML` - Open mass spectrometry data format
- Contains: Raw MS/MS spectra with m/z and intensity values

### File Matching

SpecViz automatically matches idXML files to mzML files based on:
1. **Spectrum references** in idXML (e.g., `controllerType=0 controllerNumber=1 scan=1234`)
2. **Run names** extracted from file metadata
3. **Fuzzy matching** if exact names don't match

**Example matching:**
```
sample1_ids.idXML    →  sample1_spectra.mzML
HNSCC_01.idXML       →  HNSCC_01_DDA.mzML
```

---

## idXML Processing

After uploading idXML files, SpecViz extracts and processes peptide identification data.

### What Gets Extracted

**Peptidoforms:**
- Base sequence (e.g., `PEPTIDE`)
- Modifications (e.g., `Oxidation (M)`, `Phospho (STY)`)
- Charge states
- Confidence scores (if available)

**Spectrum References:**
- File name (which mzML contains the spectrum)
- Scan number or spectrum ID
- Retention time (if available)

### Peptidoform Normalization

Multiple PSMs (Peptide-Spectrum Matches) may reference the same peptide with different modifications or charge states. SpecViz groups these into unique peptidoforms:

```
PEPTIDE+Oxidation(M)@3 [2+] → Peptidoform 1
PEPTIDE+Oxidation(M)@3 [3+] → Peptidoform 2
PEPTIDE                [2+] → Peptidoform 3
```

### Filtering Options

- **By modification**: Show only peptides with specific PTMs
- **By charge state**: Filter by precursor charge
- **By score**: Keep only high-confidence IDs (if scores available)
- **By file**: Select peptides from specific runs

---

## Peptide Selection

Choose which peptidoforms to analyze.

### Selection Methods

**1. Select All**
- Analyzes all extracted peptidoforms
- Good for: Comprehensive analysis, QC checks

**2. Manual Selection**
- Click checkboxes next to specific peptides
- Good for: Focused analysis, validation of specific IDs

**3. Search/Filter**
- Use text box to filter by sequence or modification
- Good for: Finding specific peptides in large lists

**4. Paste List**
- Copy/paste list of sequences from Excel or text file
- Automatically matches against extracted peptides
- Good for: Reproducible analysis, predefined lists

### What Happens Next

For each selected peptidoform:
1. Associated spectrum IDs are retrieved from idXML
2. Corresponding mzML files are identified
3. Spectra are loaded from mzML (indexed access, memory-efficient)
4. Spectrum metadata is extracted (precursor m/z, charge, RT)

---

## mzML Loading

SpecViz loads spectra on-demand using indexed access for memory efficiency.

### Loading Strategy

**Traditional approach** (memory-intensive):
```
Load entire mzML → Filter spectra → Process → Discard unused
```

**SpecViz approach** (memory-efficient):
```
Index mzML → Load only needed spectra → Process → Done
```

### Advantages

- Handles 10+ GB mzML files with <500 MB RAM
- Fast loading (indexed binary search)
- No temporary files needed
- Supports multiple mzML files simultaneously

### Spectrum Information Extracted

For each loaded spectrum:
- **m/z array**: Mass-to-charge ratios of fragment ions
- **Intensity array**: Signal intensities
- **Precursor info**: Selected ion m/z, charge, isolation window
- **Metadata**: Scan number, retention time, activation type

---

## Analysis Configuration

Configure fragmentation prediction and visualization settings.

### AI Model Selection

Choose prediction model based on your instrument:

**MS2PIP (Default)**
- Best for: timsTOF, modern instruments
- Model: `ms2pip_timsTOF2024`
- Speed: ~1s per spectrum
- Recommended for most users

**Prosit**
- Best for: Orbitrap (HCD fragmentation)
- Model: `Prosit_2019_intensity`
- Speed: ~2s per spectrum
- High accuracy for Orbitrap data

**AlphaPept**
- Best for: Generic predictions, unknown instruments
- Model: `AlphaPept_ms2_generic`
- Speed: ~2s per spectrum
- Good baseline for any instrument

### Fragment Settings

**Ion Types:**
- b-ions, y-ions (most common)
- a-ions, c-ions, x-ions, z-ions (less common)
- Precursor ions
- Neutral losses (H₂O, NH₃)

**Mass Tolerance:**
- Default: 0.02 Da (20 ppm for high-res)
- Adjust based on instrument accuracy
- Used for matching experimental ↔ predicted peaks

**Charge States:**
- Maximum fragment charge (default: 2)
- Usually ≤ precursor charge

**Collision Energy:**
- Default: 35 (NCE for Orbitrap)
- Adjust based on instrument method
- Affects prediction intensities

---

## Automated Analysis

SpecViz offers three analysis modes for comparing spectra.

### Mode 1: Exp vs Pred

**Purpose:** Compare experimental spectrum to AI prediction

**How it works:**
1. Loads experimental spectrum from mzML
2. Sends peptide sequence to Koina API
3. Receives predicted spectrum
4. Creates mirror plot (exp on top, pred on bottom)
5. Calculates matching metrics (R², matched ions)

**When to use:**
- Validating identifications
- Evaluating AI model performance
- Understanding fragmentation patterns

**Output:**
- Mirror plot with annotated peaks
- R² correlation score
- Number of matched ions
- Excel export with peak lists

### Mode 2: All-vs-All (Exp vs Exp)

**Purpose:** Compare all experimental spectra against each other

**How it works:**
1. Loads all selected spectra
2. Compares each spectrum to every other (n×n matrix)
3. Calculates pairwise correlations
4. Creates heatmap of similarities

**When to use:**
- Replicate quality control
- Identifying outliers
- Assessing reproducibility

**Output:**
- Correlation heatmap
- Pairwise comparison table
- Excel export with all metrics

### Mode 3: One-vs-All

**Purpose:** Compare one reference spectrum to all others

**How it works:**
1. Select reference spectrum
2. Compares reference to each other spectrum
3. Ranks by similarity
4. Identifies best matches

**When to use:**
- Finding best replicate
- Troubleshooting poor IDs
- Comparing across samples

**Output:**
- Ranked list with scores
- Individual mirror plots
- Excel export

---

## Visualization

SpecViz creates interactive mirror plots for spectrum comparison.

### Mirror Plot Anatomy

```
         Experimental (Top)
    ────────────────────────────
         m/z →
    ────────────────────────────
         Predicted (Bottom)
```

**Components:**
- **Peaks**: Vertical lines (height = intensity)
- **Colors**: Ion type (blue=b-ions, red=y-ions, etc.)
- **Labels**: m/z and ion annotation
- **Matched peaks**: Bold/highlighted
- **Unmatched peaks**: Gray/transparent

### Interactive Features

- **Zoom**: Click and drag to zoom region
- **Pan**: Shift+drag to pan view
- **Hover**: Show exact m/z and intensity
- **Legend**: Click to hide/show ion types
- **Export**: Save as PNG, SVG, or HTML

### Customization

**Plot Size:**
- Width: 800-1600 px
- Height: 400-800 px

**Labels:**
- Font size: 8-16 pt
- Show/hide ion annotations
- Show/hide m/z values

**Colors:**
- Predefined schemes (blue/red, rainbow)
- Custom colors per ion type

---

## Export Options

Save your analysis results in multiple formats.

### Excel Workbook

**Contents:**
- **Overview**: Summary table (peptides, scores, files)
- **Peak Lists**: Experimental and predicted peaks
- **Metrics**: R², matched ions, TIC
- **Settings**: Analysis parameters used

**When to use:**
- Detailed data analysis in Excel
- Statistical analysis in R/Python
- Archiving results

### HTML Interactive Plots

**Contents:**
- Standalone HTML file with embedded plot
- Fully interactive (zoom, pan, hover)
- No internet required
- Can be shared via email

**When to use:**
- Sharing results with collaborators
- Including in presentations
- Web-based reports

### PNG/SVG Images

**Contents:**
- Static image of plot
- High resolution (300 DPI for PNG)
- Vector graphics (SVG for publications)

**When to use:**
- Figures for papers/posters
- PowerPoint presentations
- Quick preview

---

## File Formats

### idXML Format

OpenMS XML format for peptide identifications.

**Structure:**
```xml
<IdentificationRun>
  <SearchParameters>...</SearchParameters>
  <ProteinIdentification>...</ProteinIdentification>
  <PeptideIdentification>
    <PeptideHit sequence="PEPTIDE" charge="2">
      <UserParam name="MS:1001117" value="50.5"/>
      <UserParam name="spectrum_reference" value="scan=1234"/>
    </PeptideHit>
  </PeptideIdentification>
</IdentificationRun>
```

**Key Elements:**
- `PeptideHit`: Individual PSM
- `sequence`: Amino acid sequence
- `charge`: Precursor charge state
- `spectrum_reference`: Link to mzML spectrum

### mzML Format

PSI standard for mass spectrometry data.

**Structure:**
```xml
<mzML>
  <run id="sample1">
    <spectrumList count="10000">
      <spectrum id="scan=1234" index="1233">
        <cvParam name="MS:1000511" value="2"/>
        <binaryDataArray>...</binaryDataArray>
      </spectrum>
    </spectrumList>
  </run>
</mzML>
```

**Key Elements:**
- `spectrum`: Individual MS/MS scan
- `binaryDataArray`: m/z and intensity arrays (compressed)
- `index`: Used for fast random access

---

## AI Models

Detailed information about supported fragmentation prediction models.

### MS2PIP

**Description:** Machine learning model trained on timsTOF data

**Features:**
- Fast predictions (~1s per spectrum)
- Excellent for modern high-resolution instruments
- Supports common modifications (Oxidation, Phospho, etc.)

**Best for:**
- timsTOF Pro/fleX
- Orbitrap Exploris
- Any high-resolution instrument

**Limitations:**
- May underperform on older low-res data
- Limited exotic modification support

### Prosit

**Description:** Deep learning model for HCD fragmentation

**Features:**
- High accuracy for Orbitrap HCD
- Trained on large datasets
- Collision energy-aware

**Best for:**
- Orbitrap Fusion, Lumos
- Q Exactive series
- Standard HCD experiments

**Limitations:**
- Slower than MS2PIP
- Designed specifically for HCD

### AlphaPept

**Description:** Generic fragmentation predictor

**Features:**
- Works across instrument types
- Good baseline performance
- Handles diverse modifications

**Best for:**
- Unknown instrument types
- Mixed instrument datasets
- Exploratory analysis

**Limitations:**
- Lower accuracy than specialized models
- May miss instrument-specific patterns

---

## Troubleshooting

Common issues and solutions.

### File Loading Issues

**Problem:** "No matching mzML files found"

**Causes:**
- mzML and idXML from different experiments
- File naming mismatch
- Files in different directories

**Solutions:**
- Check spectrum_reference in idXML matches mzML scan IDs
- Use server browser to verify file locations
- Rename files to match base names

**Problem:** "File too large" error

**Causes:**
- Browser upload limit exceeded (>2GB)
- Insufficient browser memory

**Solutions:**
- Use server browser instead
- Use Dataset Filter workspace to reduce file size
- Process files in smaller batches

### Analysis Issues

**Problem:** "Koina API timeout"

**Causes:**
- API server overloaded
- Network connectivity issues
- Too many simultaneous requests

**Solutions:**
- Retry after a few minutes
- Reduce number of spectra (< 100 at once)
- Switch to different AI model

**Problem:** "Poor R² scores" (< 0.5)

**Causes:**
- Wrong AI model for instrument
- Incorrect peptide sequence
- Poor quality spectrum

**Solutions:**
- Try different AI model (MS2PIP vs Prosit)
- Verify peptide ID confidence
- Check spectrum quality (TIC, # peaks)

### Performance Issues

**Problem:** "Out of memory" error

**Causes:**
- Too many spectra loaded simultaneously
- Large mzML files (>10GB)
- Browser memory limits

**Solutions:**
- Process fewer files at once
- Use Dataset Filter to extract relevant spectra
- Close other browser tabs/applications

**Problem:** "Slow loading times"

**Causes:**
- Network drive access
- Large number of files
- Unoptimized mzML (not indexed)

**Solutions:**
- Copy files to local drive
- Process in smaller batches
- Re-save mzML with proper indexing

---

## Advanced Features

### Custom Modification Support

Add custom modifications not in default list:

1. Go to Settings workspace
2. Click "Add Custom Modification"
3. Enter: Name, mass delta, residue specificity
4. Save and restart analysis

**Example:**
- Name: `TMT6plex`
- Mass: `+229.163`
- Residues: `K, N-term`

### Batch Processing

Process multiple files automatically:

1. Use server browser
2. Select multiple idXML files (Shift+Click)
3. Enable "Batch mode" checkbox
4. Configure once, applies to all files
5. Export generates one Excel per file

### Programmatic Access

SpecViz core functions can be imported:

```python
from core.spectrum_engine import compare_spectra
from core.mzml_processor import load_spectrum

# Load spectrum
exp_spec = load_spectrum("file.mzML", scan=1234)

# Compare to prediction
result = compare_spectra(exp_spec, pred_spec, tol=0.02)
print(f"R²: {result['r_squared']:.3f}")
```

---

## Technical Details

### Memory Optimization

SpecViz uses several strategies to minimize memory usage:

**Indexed mzML Access:**
- Only loads requested spectra (not entire file)
- Uses binary search on mzML index
- O(log n) access time

**Lazy Loading:**
- Data loaded on-demand
- Cached for repeated access
- Cleared when workspace reset

**Streaming Processing:**
- Processes spectra one at a time
- No large arrays in memory
- Constant memory footprint

### Performance Benchmarks

Typical performance on standard hardware (16GB RAM, 8 cores):

| Task | File Size | Time | Memory |
|------|-----------|------|--------|
| Load idXML | 10 MB | 2s | 50 MB |
| Load 100 spectra | 5 GB mzML | 10s | 200 MB |
| AI prediction (MS2PIP) | 1 spectrum | 1s | 100 MB |
| Generate 10 mirror plots | - | 5s | 150 MB |
| Export Excel | 100 spectra | 10s | 200 MB |

### API Rate Limits

Koina API has usage limits:

- **Free tier**: 100 requests/minute
- **Registered**: 1000 requests/minute
- **Rate limit exceeded**: Wait 60s and retry

SpecViz automatically:
- Batches predictions when possible
- Retries with exponential backoff
- Shows progress bar during API calls

---

## References

**SpecViz** is built on the OpenMS framework and WebApps template:

1. Müller, T. D., Siraj, A., et al. (2025). OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. *Journal of Proteome Research*. [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)

2. Röst, H. L., et al. (2016). OpenMS: a flexible open-source software platform for mass spectrometry data analysis. *Nature Methods*, 13, 741–748.

3. Gabriels, R., et al. (2019). Updated MS²PIP web server delivers fast and accurate MS² peak intensity prediction for multiple fragmentation methods. *Nucleic Acids Research*, 47(W1), W295-W299.

4. Gessulat, S., et al. (2019). Prosit: proteome-wide prediction of peptide tandem mass spectra by deep learning. *Nature Methods*, 16, 509–518.
