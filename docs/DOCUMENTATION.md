# MSpecViz Documentation

Comprehensive reference for every feature, setting, and workflow step in MSpecViz — a Streamlit-based application for MS/MS spectrum visualization and comparison.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation & Access](#installation--access)
3. [Workflow Overview](#workflow-overview)
4. [Step 1: File Selection](#step-1-file-selection)
5. [Step 2: Target Peptides](#step-2-target-peptides)
6. [Step 3: Reference Identification](#step-3-reference-identification)
7. [Step 4: Settings](#step-4-settings)
8. [Step 5: Search Peptides](#step-5-search-peptides)
9. [Step 6: Process & Export](#step-6-process--export)
10. [Step 7: Results](#step-7-results)
11. [Parquet Fast Loading](#parquet-fast-loading)
12. [Caching System](#caching-system)
13. [Comparison Metrics](#comparison-metrics)
14. [File Formats](#file-formats)
15. [Configuration](#configuration)
16. [Troubleshooting](#troubleshooting)
17. [Architecture & Technical Details](#architecture--technical-details)

---

## Overview

MSpecViz compares experimental MS/MS spectra against reference spectra for target peptides. It annotates peaks with theoretical fragment ions and generates publication-ready mirror plots showing how well two spectra agree.

### What MSpecViz Does

1. **Searches** idXML identification files for peptide-spectrum matches (PSMs) of your target peptides
2. **Loads** the corresponding MS/MS spectra from mzML files
3. **Annotates** experimental peaks by matching them to theoretical fragment ions (b, y, a, c, x, z, precursor, immonium) within a configurable m/z tolerance
4. **Compares** query spectra against reference spectra using seven quantitative metrics
5. **Generates** mirror plots (query on top, reference on bottom) with color-coded ion annotations
6. **Exports** plots as PNG/SVG images and comparison metrics as CSV/XLSX spreadsheets

### Typical Use Cases

- **Validating peptide identifications** — Compare experimental fragmentation against a reference standard (e.g., synthetic peptide library)
- **Quality control** — Check replicate consistency by comparing spectra from different runs
- **Method development** — Evaluate how well a fragmentation method reproduces expected patterns
- **Publication figures** — Generate annotated mirror plots for manuscripts and presentations

---

## Installation & Access

### Local Installation

```bash
# Create the conda environment (first time only)
conda env create -f environment.yml
conda activate mspecviz

# Start the app
streamlit run app.py
```

The app opens at http://localhost:8501.

### Remote Access (SSH Tunnel)

If MSpecViz runs on a remote server accessible only via SSH through a jumphost:

```bash
ssh -L 8501:localhost:8501 -J jumphost server_name
```

Then open http://localhost:8501 in your local browser.

A convenience script is provided:

```bash
./scripts/connect_mspecviz.sh                    # defaults
./scripts/connect_mspecviz.sh myuser@jumphost myuser@server  # custom
LOCAL_PORT=9501 ./scripts/connect_mspecviz.sh     # different local port
```

---

## Workflow Overview

MSpecViz is organized into **7 sequential tabs**. You work through them left to right:

```
1. Files  →  2. Peptides  →  3. Reference  →  4. Settings  →  5. Search  →  6. Process  →  7. Results
```

Each step validates its prerequisites before allowing you to proceed. You can always go back to earlier tabs to change settings.

**Shortcut**: If you have a pre-converted Parquet file (from a previous run or the CLI tool), select it in Step 1 and skip directly to Step 6.

---

## Step 1: File Selection

Select the input files for your analysis. See [File Selection](file_selection.md) for full details.

### Directory Browser

Each file type (idXML and mzML) has its own directory browser with:

- **Path input**: Type or paste a directory path directly
- **Up button**: Navigate to the parent directory
- **Folder grid**: Click folders to navigate into them (up to 15 shown)
- **File filter**: Type to filter files by name
- **Select All**: Checkbox to select all matching files in the current directory
- **Multi-select**: Pick individual files from the list

### Accumulating Selections

You can select files from multiple directories. Selections persist as you navigate — files from previous directories remain selected. The "All selected files" expander shows everything grouped by directory with individual remove buttons.

### Same-Name Files

If you select files with identical names from different directories, MSpecViz automatically disambiguates them by prepending the minimum number of parent directory components needed to make each name unique. For example:

```
/data/experiment1/sample.mzML  →  experiment1/sample.mzML
/data/experiment2/sample.mzML  →  experiment2/sample.mzML
/data/unique_file.mzML         →  unique_file.mzML
```

Exact duplicate paths (same file selected twice) are silently removed.

### Parquet Shortcut

Below the file browsers, you can paste the path to a pre-converted `.parquet` file. If valid, this skips Steps 2–5 entirely — the Parquet file already contains all the spectra, peptide information, and reference categorization.

### Validation

Step 1 requires either:
- At least 1 idXML **and** 1 mzML file selected, **or**
- A valid Parquet file path

---

## Step 2: Target Peptides

Enter the peptide sequences you want to analyze. See [Peptide Search](peptide_search.md) for details on how the search works.

### Input Format

Enter one peptide per line in the text area:

```
LGGNEQVTR
GAGSSEPVTGLDAK
VEATFGVDESNAK
```

**Important**: Enter **base sequences** (without modifications). MSpecViz automatically finds all modified variants (e.g., Oxidation, Phosphorylation) during the search.

- Case insensitive — `peptidek` and `PEPTIDEK` are treated the same
- Whitespace and empty lines are ignored
- The counter below the text area confirms how many peptides you've entered

### Preview

Click **Preview in idXML** to scan your selected idXML files and display a sample of peptides they contain. This helps verify that your target peptides exist in the data before running the full search.

---

## Step 3: Reference Identification

Configure how MSpecViz separates your mzML files into **reference** and **query** groups. See [Reference Setup](reference_setup.md) for detailed explanation.

### How It Works

Enter a **substring** that identifies your reference files. MSpecViz checks each mzML filename:

- If the filename **contains** the substring (case-insensitive) → **Reference**
- Otherwise → **Query**

**Example**: If your reference files are named `Reference_01.mzML`, `Reference_02.mzML` and your query files are `Sample_A.mzML`, `Sample_B.mzML`, enter `Reference`.

### Preview

After entering a substring, a two-column preview shows which files are categorized as Reference and which as Query. Verify this looks correct before proceeding.

### No Reference (Optional)

If you leave the reference substring blank, MSpecViz generates **individual annotated spectrum plots** instead of mirror plots. Each spectrum is plotted on its own with theoretical fragment annotations but no comparison. This is useful for simply visualizing fragmentation patterns.

### Validation

If a substring is entered, at least one file must match as Reference and at least one must remain as Query.

---

## Step 4: Settings

Configure the annotation and visualization parameters. See [Annotation Settings](annotation_settings.md) for detailed explanations.

### Ion Types

Select which fragment ion series to annotate. Available types:

| Ion Type | Description | Default |
|----------|-------------|---------|
| **b** | N-terminal fragment ions (CID/HCD) | Yes |
| **y** | C-terminal fragment ions (CID/HCD) | Yes |
| **a** | N-terminal minus CO | No |
| **c** | N-terminal plus NH₃ (ETD) | No |
| **x** | C-terminal plus CO (ETD) | No |
| **z** | C-terminal minus NH₃ (ETD) | No |
| **precursor** | Intact precursor ion at various charges | No |
| **immonium** | Single amino acid internal fragments | No |

### m/z Tolerance

How close (in Daltons) an experimental peak must be to a theoretical fragment m/z to be considered a match. Default: **0.02 Da**. Increase for lower-resolution data, decrease for ultra-high-resolution instruments.

### Annotation Options

| Option | Default | What It Does |
|--------|---------|--------------|
| **Neutral losses** | Off | Include H₂O and NH₃ loss peaks in the theoretical spectrum |
| **First prefix ion** | On | Include the first prefix ion (a₁ or b₁), which some search engines exclude |
| **Remove precursor peaks** | Off | Remove all precursor ions from the plot. Intensities are renormalized afterward |

### Intensity Mode

- **Relative** (default): Each spectrum is normalized to its base peak (0–1 scale). Best for comparing fragmentation patterns regardless of signal intensity.
- **Raw**: Original intensity values from the mzML file. Shows actual signal strength but makes cross-spectrum comparison harder.

---

## Step 5: Search Peptides

Click **Search Peptides** to find all PSMs matching your target peptides. See [Peptide Search](peptide_search.md) for the full algorithm.

### What Happens During Search

1. Each idXML file is loaded with PyOpenMS
2. All `PeptideIdentification` records are scanned
3. PSMs whose base sequence matches any target peptide are collected
4. The `SpectraData` metadata in the idXML maps each PSM to a specific mzML file
5. If multiple mzML files share the same name, MSpecViz uses path suffix matching and directory proximity heuristics to resolve ambiguity

### Results Table

The results table shows:

| Column | Meaning |
|--------|---------|
| **Peptide** | Base amino acid sequence |
| **Modification** | Post-translational modifications found (e.g., Oxidation@M5) |
| **PSMs** | Number of peptide-spectrum matches |
| **idXML Files** | Number of idXML files containing this peptide |
| **mzML Resolved** | Number of PSMs successfully mapped to mzML files |
| **Unresolved** | PSMs that could not be mapped (missing metadata) |

### Convert to Parquet (Optional)

After searching, you can convert the matched spectra to a Parquet file for instant loading on future runs. Enter an output path and click **Convert to Parquet**. See [Parquet Fast Loading](#parquet-fast-loading).

---

## Step 6: Process & Export

Load spectra, generate plots, and export everything. See [Export](export.md) for output format details and [Mirror Plots](mirror_plots.md) for plot interpretation.

### Export Settings

| Setting | Default | Description |
|---------|---------|-------------|
| **Project Name** | `MSpecViz_analysis` | Base name for the export folder |
| **Export Directory** | `./exports` | Where export folders are created |
| **Export Formats** | PNG | Image formats to generate (PNG and/or SVG) |
| **Plot unpaired** | Off | Whether to plot spectra that have no matching reference |

### Processing Pipeline

When you click **Process & Export**, MSpecViz executes the following pipeline:

1. **Cache check**: Computes a fingerprint of your inputs and checks for an existing Parquet cache
2. **Load references**: All reference spectra for each peptide variant are loaded first and cached in memory
3. **Load queries**: Query spectra are loaded one at a time (streaming)
4. **Annotate**: Each spectrum's peaks are matched against theoretical fragments using PyOpenMS's `SpectrumAlignment`
5. **Parallel plot generation**: All annotated spectrum pairs are processed in parallel using multiple worker processes. Each worker independently computes metrics, creates a mirror plot (or single-spectrum plot), and saves it to disk. The number of workers is configurable (default: 10 in Streamlit, 4 in CLI).
6. **Summary export**: A comparison summary spreadsheet (CSV + XLSX) is written with per-pair metrics
7. **Cache save**: If this was a fresh mzML load, a Parquet cache is saved for next time

### Progress Tracking

Two separate progress bars track:
- **Reference spectra**: Loading and counting references
- **Query spectra**: Processing, plotting, and exporting each query

Both show elapsed time and estimated remaining time.

### Unpaired Spectra

When a query spectrum has no corresponding reference (e.g., the peptide was found only in query files), the behavior depends on the **Plot unpaired** checkbox:

- **Off (default)**: The spectrum is silently skipped. Only mirror plots are generated.
- **On**: A single annotated spectrum plot is created (no mirror, just the query with annotations).

---

## Step 7: Results

Displays a summary of the completed analysis:

- **Counts**: Number of reference spectra, query spectra, and plots generated
- **Files**: PNG and SVG file counts
- **Summary files**: Links to `comparison_summary.csv` and `comparison_summary.xlsx`
- **Export location**: Full path to the export directory
- **Errors**: Any issues encountered during processing (expandable section)

Expandable sections list individual exported files by format.

---

## Parquet Fast Loading

MSpecViz supports loading spectra from pre-converted Parquet files instead of parsing mzML XML. This is 100–1000x faster.

### Three Ways to Use Parquet

1. **Auto-cache**: After a successful run, MSpecViz automatically saves a Parquet cache. Next time you use the same inputs, it loads from cache instead of mzML.

2. **In-app conversion**: After Step 5 (Search), click "Convert to Parquet" to save matched spectra. Use the Parquet file path in Step 1 on future runs.

3. **CLI conversion**: For large datasets, use the standalone conversion tool:
   ```bash
   conda run -n mspecviz python tools/convert_mzml_to_parquet.py \
     --idxml /path/to/*.idXML \
     --mzml-dir /path/to/mzml_files/ \
     --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
     --output exports/my_spectra.parquet \
     --reference-substring "Synpep"
   ```
   See [tools/README.md](../tools/README.md) for full CLI documentation.

### What's Stored in a Parquet File

Each row represents one spectrum:

| Column | Type | Description |
|--------|------|-------------|
| `source_mzml` | string | Source mzML display name |
| `spectrum_id` | string | Native spectrum ID (e.g., `index=12345`) |
| `peptidoform` | string | Full sequence with modifications |
| `base_sequence` | string | Unmodified peptide sequence |
| `variant_key` | string | Grouping key (`sequence\|modifications`) |
| `charge` | int | Precursor charge state |
| `precursor_mz` | float | Precursor m/z |
| `rt` | float | Retention time |
| `is_reference` | bool | Whether this is a reference spectrum |
| `mz` | float[] | m/z array |
| `intensity` | float[] | Intensity array |
| `score` | float | PSM score |
| `idxml_file` | string | Source idXML display name |

### Performance Comparison

| Scenario | mzML (first run) | mzML (cached index) | Parquet |
|----------|-------------------|---------------------|---------|
| 10 spectra, 1 file (100k spectra) | 5–30s | <1s | ~0.1s |
| 100 spectra, 7 files | 35–210s | ~4s | ~0.5s |
| Re-run same analysis | Same as first | Same as first | ~0.5s |

---

## Caching System

MSpecViz automatically caches loaded spectra in Parquet format so that repeated runs with the same inputs load instantly. See [Caching](caching.md) for the full technical reference.

### Key Points

- A **fingerprint** is computed from your inputs (file paths + modification times + sizes, peptide list, reference substring)
- If the fingerprint matches an existing cache file, spectra are loaded from Parquet instead of mzML
- Cache is saved **only after a complete, successful run** — interrupted runs leave no cache
- Cache is automatically invalidated when any input changes
- After code changes that affect spectrum loading, delete the cache manually:
  ```bash
  rm exports/.mspecviz_cache/*.parquet exports/.mspecviz_cache/*.json
  ```

---

## Comparison Metrics

Seven quantitative metrics are computed for each query-vs-reference pair. See [Spectrum Metrics](spectrum_metrics.md) for mathematical details, equations, and interpretation guides.

### Quick Reference

| Metric | Abbreviation | What It Answers |
|--------|-------------|-----------------|
| R² Matching | R2m | Do shared annotated ions have similar relative intensities? |
| Cosine Similarity | Cos | How similar are the intensity vectors of shared ions? |
| Dice Coefficient | Dice | What fraction of annotated ions are shared? |
| Jaccard Index | Jacc | Same as Dice but with stricter (union) denominator |
| Shared Peak Ratio | ShR | Is the smaller spectrum's information contained in the larger? |
| R² Binned | R2b | Do the overall spectral profiles correlate? (all peaks) |
| Spectral Entropy | Ent | Information-theoretic similarity (all peaks) |

### Output

Metrics appear in two places:
1. **On each mirror plot**: Abbreviated metrics text in the top-right corner
2. **Summary spreadsheet**: `comparison_summary.csv` and `.xlsx` in the export directory, with one row per query-reference pair

---

## File Formats

### Input: idXML

OpenMS XML format for peptide identifications. Contains:

- **PeptideIdentification** records: Each contains one or more `PeptideHit` with sequence, modifications, charge, and score
- **SpectraData** metadata: Maps spectrum references to mzML files via `id_merge_index`
- Created by OpenMS search engines (MSGF+, Comet, etc.) or converters

### Input: mzML

PSI standard format for mass spectrometry data. Contains:

- **Spectrum records**: MS/MS scans with m/z and intensity arrays (binary-compressed)
- **Index**: Enables fast random access to individual spectra without loading the entire file
- **Native IDs**: Spectrum identifiers (e.g., `index=12345`, `scan=6789`)

**Note on indexing**: Some instruments (Bruker/timsTOF) use 1-based native IDs (`index=1` is the first spectrum). MSpecViz auto-detects the offset from the first spectrum in each file.

### Output: PNG/SVG

- **PNG**: Raster images at 300 DPI. Good for presentations, reports, quick review.
- **SVG**: Vector graphics. Good for publications (infinitely scalable, editable in Illustrator/Inkscape).

### Output: CSV/XLSX

The comparison summary contains one row per query-reference pair with:
- Shortened filenames for readability
- Peptide sequence and modification
- Charge state
- All 7 comparison metrics
- Peak count context (total peaks, annotated peaks, shared ions)
- Full filenames

---

## Configuration

All settings are in `settings.json`. Changes take effect on the next app restart.

### Key Settings

```json
{
  "data_directories": {
    "default": "/"                    // Starting directory for file browser
  },
  "ion_types": {
    "default": ["b", "y"],            // Default ion types selected
    "available": ["a","b","c","x","y","z","precursor","immonium"]
  },
  "mz_tolerance": {
    "default": 0.02,                  // Default m/z tolerance in Da
    "unit": "Da"
  },
  "annotation": {
    "add_losses": false,              // Neutral losses (H₂O, NH₃)
    "add_first_prefix_ion": true,     // Include a₁/b₁
    "remove_precursor_peaks": false   // Strip precursor ions
  },
  "export": {
    "default_directory": "./exports", // Export base directory
    "project_name": "MSpecViz_analysis",
    "formats": ["png"],              // Default export formats
    "dpi": 300,                      // PNG resolution
    "figsize": [12, 8]               // Plot size (width, height inches)
  },
  "performance": {
    "workers": 4                     // Parallel worker processes for plot generation
  },
  "cache": {
    "enabled": true,                 // Auto-caching on/off
    "subdirectory": ".mspecviz_cache" // Cache folder name
  },
  "ion_colors": {                    // Hex colors for each ion type
    "b": "#0066CC",
    "y": "#CC0000",
    ...
  }
}
```

See [Annotation Settings](annotation_settings.md) for detailed explanations of each annotation option.

---

## Troubleshooting

### File Selection Issues

**No files appear in the directory browser**
- Verify the directory exists and is readable
- Check that files have the expected extension (`.idXML`, `.mzML` — case-sensitive on Linux)
- Try typing the full path in the path input box

**Files from different directories have the same name**
- MSpecViz handles this automatically by adding parent directory components to disambiguate
- Exact duplicate paths (same file selected twice) are removed silently

### Search Issues

**"No matches found" for peptides that should be there**
- Use "Preview in idXML" to check which peptides exist in your files
- Verify you're entering base sequences without modifications
- Check that idXML files are from the expected experiment

**PSMs could not be mapped to mzML files**
- The idXML's `SpectraData` metadata may be missing or point to mzML files not in your selection
- Ensure you've selected all mzML files referenced by the idXML
- Check the `id_merge_index` entries in the idXML's `ProteinIdentification` section

### Processing Issues

**Annotation errors**
- Empty m/z arrays: The spectrum may be corrupt or empty in the mzML file
- Invalid charge: The PSM may have charge=0 (missing charge information)
- Parse errors: The peptidoform string may contain unsupported modification notation

**Slow first run**
- Normal for large mzML files. MSpecViz builds an index and streams spectra.
- Subsequent runs with the same inputs use the Parquet cache (~100x faster)
- For very large datasets, pre-convert to Parquet using the CLI tool

**Cache shows wrong results after a code update**
- The fingerprint only tracks input files, not code. Delete the cache:
  ```bash
  rm exports/.mspecviz_cache/*.parquet exports/.mspecviz_cache/*.json
  ```

### Connection Issues (Remote Server)

**"localhost refused to connect"**
- You need an SSH tunnel. See [Remote Access](#installation--access)
- Verify Streamlit is actually running on the server: `curl -s -o /dev/null -w '%{http_code}' http://localhost:8501` (should return 200)
- Check that port 8501 isn't already in use locally

---

## Architecture & Technical Details

### Processing Pipeline

```
idXML files ──→ search_peptides_in_idxmls() ──→ VariantGroups
                                                      │
mzML files  ──→ load_spectra_streaming() ─────────────┤
                      │                                │
                      ├── open_indexed_mzml()          │
                      ├── build_mzml_to_matches_map()  │
                      └── stream spectra by variant ───┘
                                                       │
                              ┌─────────────────────────┘
                              │
                    prepare_spectrum_for_plot()
                              │
                    ├── generate_theoretical_fragments()
                    │         (PyOpenMS TheoreticalSpectrumGenerator)
                    │
                    ├── annotate_spectrum()
                    │         (PyOpenMS SpectrumAlignment)
                    │
                    └── normalize / filter
                              │
                    ──→  ComparisonTask objects (collected)
                              │
                    run_parallel_comparisons()
                              │
                    ├── create_mirror_plot()  ──→  PNG/SVG export
                    └── compute_all_metrics()  ──→  CSV/XLSX summary
```

### Key Backend Modules

| Module | Responsibility |
|--------|---------------|
| `src/peptide_search.py` | idXML parsing, PSM extraction, mzML resolution |
| `src/spectrum_loader.py` | mzML indexed access, Parquet loading, streaming |
| `src/peak_annotation.py` | Theoretical fragment generation, peak matching |
| `src/mirror_plot.py` | Matplotlib mirror plot and single spectrum plots |
| `src/spectrum_metrics.py` | Seven comparison metrics |
| `src/parallel_processor.py` | Parallel plot generation with ProcessPoolExecutor |
| `src/parquet_cache.py` | Fingerprinting, cache lookup, cache save |
| `src/path_utils.py` | Path disambiguation, deduplication |
| `src/export.py` | Filename sanitization, directory creation |

### Memory Strategy

MSpecViz uses **streaming processing** to minimize memory usage:

- mzML files are opened with `OnDiscMSExperiment` (indexed access, no full load)
- Reference spectra are loaded per-variant and held in memory
- Query spectra are loaded one at a time, annotated, and collected as lightweight task objects
- Plot generation runs in parallel worker processes; each worker creates one figure, saves it, and releases memory immediately
- This enables processing large datasets (10+ GB mzML files) with modest RAM

### mzML Index Handling

MSpecViz supports multiple native ID formats:
- `index=N` (common in Bruker/timsTOF data)
- `scan=N` (common in Thermo data)
- `spectrum=N`

For `index=` IDs, the offset between native ID and array position is auto-detected from the first spectrum in each file. Bruker/timsTOF files typically use 1-based indexing (`index=1` is the first spectrum at array position 0).

---

## References

1. Muller, T. D., Siraj, A., et al. (2025). OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. *Journal of Proteome Research*. [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)

2. Rost, H. L., et al. (2016). OpenMS: a flexible open-source software platform for mass spectrometry data analysis. *Nature Methods*, 13, 741-748.
