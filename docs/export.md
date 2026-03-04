# Export & Output

This document describes everything MSpecViz exports: plots, summary spreadsheets, and the directory structure.

---

## Export Directory Structure

Each run creates a timestamped directory under your chosen export path:

```
exports/
└── MSpecViz_analysis_20260219_143000/
    ├── comparison_summary.csv          # Metrics for every query-reference pair
    ├── comparison_summary.xlsx         # Same data in Excel format
    ├── png/                            # PNG plot images
    │   ├── LGGNEQVTR/                 # One folder per peptide variant
    │   │   ├── LGGNEQVTR_SampleA_vs_Reference01.png
    │   │   ├── LGGNEQVTR_SampleA_vs_Reference02.png
    │   │   └── LGGNEQVTR_SampleB_vs_Reference01.png
    │   └── PEPTMIDEK_M5_Ox_/
    │       └── PEPTMIDEK_SampleA_vs_Reference01.png
    └── svg/                            # SVG plots (if enabled)
        ├── LGGNEQVTR/
        │   └── ...
        └── PEPTMIDEK_M5_Ox_/
            └── ...
```

### Folder Naming

The root folder name combines:
- **Project name**: Set in Export Settings (default: `MSpecViz_analysis`)
- **Timestamp**: `YYYYMMDD_HHMMSS` format

Characters invalid for filenames are replaced with underscores. The name is truncated to 100 characters.

### Variant Subfolders

Plots are grouped by peptide variant:
- Unmodified: `PEPTIDEK/`
- Modified: `PEPTIDEK_M5_Ox_/` (compact modification key sanitized for filesystem)

---

## Plot Files

### PNG (Default)

- Resolution: 300 DPI
- Size: ~12 x 8 inches (configurable in `settings.json`)
- Background: White
- Typical file size: 200–500 KB per plot

### SVG (Optional)

- Vector format — infinitely scalable
- Editable in Illustrator, Inkscape, or any SVG editor
- Background: White
- Typical file size: 50–200 KB per plot

### Filename Convention

**Mirror plots** (query vs reference):
```
{peptide}_{query_filename_stem}_vs_{ref_filename_stem}.{ext}
```

**Single spectrum plots** (no reference):
```
{peptide}_{source_filename_stem}.{ext}
```

Filenames are sanitized:
- Invalid characters (`/ \ : * ? " < > |`) replaced with `_`
- Source filename stems truncated to 30 characters
- Total filename truncated to 100 characters

---

## Comparison Summary Spreadsheet

After processing, MSpecViz writes a comparison summary with one row per query-reference pair.

### CSV (`comparison_summary.csv`)

Always written. Can be opened in any spreadsheet application or loaded in Python/R for further analysis.

### XLSX (`comparison_summary.xlsx`)

Written if the `openpyxl` package is available. Same data as CSV but in native Excel format.

### Columns

| Column | Description |
|--------|-------------|
| **Query** | Shortened query filename (common prefixes/suffixes removed) |
| **Reference** | Shortened reference filename |
| **Peptide** | Base amino acid sequence |
| **Modification** | Modification string (empty if unmodified) |
| **Charge** | Precursor charge state of the query spectrum |
| **R² match** | R² from linear regression on shared annotated ion intensities |
| **R² binned** | R² from binned intensity profiles (all peaks) |
| **Cosine** | Cosine similarity of shared annotated ion intensities |
| **Entropy** | Spectral entropy similarity (all peaks, information-theoretic) |
| **Dice** | Dice coefficient of annotated ion overlap |
| **Jaccard** | Jaccard index of annotated ion overlap |
| **SharedRatio** | Fraction of the smaller annotated set that is shared |
| **Peaks (Q)** | Total peaks in the query spectrum |
| **Peaks (R)** | Total peaks in the reference spectrum |
| **Annotated (Q)** | Annotated peaks in the query |
| **Annotated (R)** | Annotated peaks in the reference |
| **Shared Ions** | Number of annotation labels shared between both spectra |
| **Query File** | Full query filename |
| **Ref File** | Full reference filename |

### Shortened Filenames

The "Query" and "Reference" columns use shortened names for readability. The shortening algorithm:

1. Strips file extensions (`.mzML`)
2. Removes common prefix tokens shared by all filenames
3. Removes common suffix tokens
4. Keeps the distinguishing parts

Example:
```
Exp01_Sample_A.mzML   →  Sample_A
Exp01_Sample_B.mzML   →  Sample_B
Exp01_Synpep_01.mzML  →  Synpep_01
```

The full filenames are preserved in the "Query File" and "Ref File" columns.

---

## Export Settings

Configured in Step 6 (Process & Export):

| Setting | Default | Description |
|---------|---------|-------------|
| **Project Name** | `MSpecViz_analysis` | Base name for the export folder |
| **Export Directory** | `./exports` | Parent directory for export folders |
| **Export Formats** | PNG | Which image formats to generate |
| **Plot unpaired** | Off | Export single-spectrum plots for unmatched queries |

### Format Selection

- **PNG only** (default): Fastest export, smallest disk usage
- **SVG only**: Vector output only
- **Both PNG and SVG**: Creates parallel directory trees

### Defaults in `settings.json`

```json
{
  "export": {
    "default_directory": "./exports",
    "project_name": "MSpecViz_analysis",
    "formats": ["png"],
    "dpi": 300,
    "figsize": [12, 8]
  }
}
```

---

## Understanding the Output

### How Many Plots to Expect

The number of plots depends on your data:

| Scenario | Plots Generated |
|----------|----------------|
| 5 queries x 2 references, 1 peptide | 10 mirror plots |
| 5 queries x 2 references, 3 peptides (each in all files) | 30 mirror plots |
| 10 queries, 0 references (no substring), unpaired=on | 10 single plots |
| 10 queries, 0 references, unpaired=off | 0 plots |

### Metrics Values

See [Spectrum Metrics](spectrum_metrics.md) for how to interpret metric values. A quick guide:

| Metric Range | Quality |
|-------------|---------|
| R²/Cos > 0.9 | Excellent match |
| R²/Cos 0.7–0.9 | Good match |
| R²/Cos 0.3–0.7 | Moderate match — review the plot |
| R²/Cos < 0.3 | Poor match — likely different spectra |
| Dice > 0.7 | Most annotated ions are shared |
| Dice < 0.3 | Few shared ions |

### NaN Values

A metric may be `NaN` when there isn't enough data to compute it:
- **R²**: Requires at least 2 shared ions (or 2 shared bins)
- **Cosine**: Requires at least 1 shared ion
- **SharedRatio**: Requires at least 1 annotated peak in both spectra

---

## Parquet Cache (Automatic)

In addition to the export directory, MSpecViz may save a Parquet cache for instant reloading. This is stored separately:

```
exports/.mspecviz_cache/
├── 48020b32ce9d.parquet    # Cached spectra
└── 48020b32ce9d.json       # Sidecar metadata
```

This is not part of the analysis output — it's an internal optimization. See [Caching](caching.md) for details.
