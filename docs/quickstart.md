# MSpecViz Quick Start Guide

This guide walks through the 4-tab Streamlit workflow using test data, and
shows the equivalent CLI invocation. For deeper references see:

- [architecture/overview.md](architecture/overview.md) — how the codebase is organised.
- [settings/overview.md](settings/overview.md) — every configuration key explained.
- [cli/usage.md](cli/usage.md) — full CLI surface.
- [metrics/comparison.md](metrics/comparison.md) — formulas and ranges for the
  7 comparison metrics (R2m, R2b, Cos, Ent, Dice, Jacc, ShR).
- [export/layout.md](export/layout.md) — output directory layout and file formats.

---

## What You Need

| File Type | Extension | Contains |
|---|---|---|
| Identification files | `.idXML` | Peptide-spectrum matches from a search engine |
| Spectrum files | `.mzML` | Raw MS/MS spectra (m/z and intensity arrays) |

Both file types must come from the same experiment. The idXML files reference
spectra inside the mzML files by spectrum ID.

---

## Test Data

Use these paths for a quick verification run:

| Input | Path |
|---|---|
| idXML | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/ProT001_1_pout_filtered.idXML` |
| mzML directory | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/` |
| Peptide list | `/mnt/volume/nacho/data/peptide_lists/peptides.txt` |

Expected results: 2 PSMs (ASDDFDQASY and RTEFNLNQY), 373 total peaks.

---

## Starting the App

```bash
conda activate mspecviz
cd streamlit_app
streamlit run app.py
```

Run from inside `streamlit_app/` so Streamlit picks up
`.streamlit/config.toml`. Workspaces still resolve to `<project-root>/workspaces/`,
so `cd` does not strand your existing data.

Open http://localhost:8501 in your browser.

To run headlessly without the Streamlit UI, use the CLI from the project root:

```bash
python -m cli cli/config/config_example.yaml
```

The CLI is YAML-only (no command-line flag overrides) — copy the example
and edit it for per-run tweaks. Set `MSPECVIZ_LOG_LEVEL=ERROR` for a
quiet run.

---

## Tab 1: Setup

The Setup tab configures all inputs and parameters before running the analysis.

### File Selection

1. **idXML files** -- Use the directory browser in the left column to navigate
   to your idXML files. Click file names to select them. Selections accumulate
   across directories.

2. **mzML files** -- Use the directory browser in the right column to select
   the corresponding mzML files. These must be the files referenced in the
   idXML metadata.

### Peptide Input

Enter target peptide sequences in the text area, one per line:

```
ASDDFDQASY
RTEFNLNQY
```

Enter **base sequences** (unmodified). MSpecViz automatically finds all
modified variants during the search step. Sequences are case-insensitive.

Use "Preview peptides in idXML" to see which peptides are present in your
selected idXML files before committing to a search.

### Reference Configuration

Enter a **reference substring** that identifies your reference mzML files
(for example: `Synpep`, `Standard`, `Ref`). Files whose names contain this
substring are categorized as references; all others are queries.

- With a reference substring set, MSpecViz generates mirror plots (query on
  top, reference on bottom) and computes comparison metrics.
- Without a reference substring, each spectrum is plotted individually with
  annotations only (no mirror comparison).

### Koina AI Prediction

When no experimental reference is available, enable Koina AI prediction as a
fallback. Configure:

- **Model** -- Choose from 14 models across MS2PIP, Prosit, AlphaPept, UniSpec,
  and PredFull families.
- **Collision energy** -- Normalized collision energy (default: 35).
- **Instrument type** -- QE, LUMOS, TIMSTOF, or SCIEXTOF.

Koina predicts both m/z positions and relative intensities, producing a
synthetic reference spectrum for mirror plot comparison.

### Annotation Settings

- **Ion types** -- Fragment ion series to annotate (default: b and y). Options
  include a, b, c, x, y, z, precursor, and immonium ions.
- **m/z tolerance** -- Maximum distance in Da between an experimental peak and
  a theoretical fragment for a match (default: 0.02 Da).
- **Neutral losses** -- Include H2O and NH3 loss peaks in annotation.
- **Intensity mode** -- Relative (normalized to base peak = 100%) or Raw.

### Performance Settings

- **Workers** -- Number of parallel processes for export (default: 4). Set to
  1 for sequential processing (useful for debugging).

---

## Tab 2: Search & Load

This tab runs the analysis pipeline in two stages.

### Stage 1: Search Peptides

Click **Search Peptides** to scan selected idXML files for target peptide
sequences. MSpecViz finds all PSMs, groups them by peptide variant (base
sequence plus modification state), and resolves which mzML file contains each
spectrum.

The results table shows:

| Column | Meaning |
|---|---|
| Peptide | Base sequence |
| Modification | Modification state (or "unmodified") |
| PSMs | Number of peptide-spectrum matches found |
| idXML Files | Source identification files |
| mzML Resolved | Whether the corresponding mzML files were located |

### Stage 2: Search & Load Spectra

After reviewing search results, click **Search & Load Spectra** to:

1. Stream spectra from mzML files.
2. Annotate each spectrum with theoretical fragment ions.
3. Pair query spectra with references (by filename substring match).
4. If no experimental reference exists and Koina is enabled, generate an
   AI-predicted reference spectrum.
5. Compute the 7 comparison metrics for each query/reference pair.
6. Build Polars DataFrames for the openms_insight interactive explorer.

A progress bar tracks loading. When complete, proceed to the Explore tab.

---

## Tab 3: Explore

The Explore tab offers three sub-views, selectable via radio buttons.

### Interactive Explorer

Renders four linked openms_insight components:

- **PSM Table** -- Lists all loaded PSMs with peptide, modification, charge,
  score, and source file. Click any row to select it.
- **SequenceView** -- Shows the peptide sequence with colored fragment ion
  coverage for the selected PSM.
- **LinePlot** -- Displays the annotated stick spectrum for the selected PSM.
  Peaks are colored by ion type (matching SequenceView colors). Hover for
  annotation labels and m/z values.
- **Heatmap** -- RT-by-m/z intensity map across all loaded spectra. Useful for
  spotting chromatographic patterns and confirming retention times.

All four components are linked through a shared StateManager. Selecting a row
in the Table updates the LinePlot and SequenceView automatically.

### Mirror Plots

Displays Plotly interactive mirror plots for each query-reference pair:

- **Upper panel** -- Query (experimental) spectrum, peaks pointing upward.
- **Lower panel** -- Reference spectrum, peaks pointing downward.
- Annotated peaks are colored by ion type; unannotated peaks are gray.
- Hover over any peak to see its m/z, intensity, and annotation label.
- Comparison metrics are displayed alongside the plot.

Use the PSM selector to browse different peptide variants and source files.
When no reference is available, a single-spectrum plot is shown instead.

### Comparison Metrics

Displays the full comparison data in tabular form:

- **Metrics table** -- All 7 comparison metrics for every query-reference pair.
- **Summary statistics** -- Aggregate statistics (mean, median, std) across
  all comparisons.

See [metrics/comparison.md](metrics/comparison.md) for detailed definitions
and [metrics/interpretation.md](metrics/interpretation.md) for guidance on
which metric to trust when.

---

## Tab 4: Export

The Export tab generates publication-quality outputs and manages export history.

### Export Options

- **Project name** -- Names the versioned export directory.
- **Export directory** -- Where to save outputs.
- **Formats** -- PNG (default), SVG, PDF for static plots; HTML for interactive
  Plotly plots.
- **DPI** -- Resolution for raster images (default: 300).

### What Gets Exported

- Parquet DataFrames (psm.parquet, peaks.parquet, seq.parquet, heatmap.parquet)
- CSV and XLSX comparison summary with the 7 metrics plus per-pair context
  columns (peak counts, annotation counts, shared peaks)
- Mirror plots organized by peptide variant in format-specific subdirectories
- Interactive HTML plots (Plotly, using CDN for small file sizes)
- Summary analysis plots (dual-axis bar/box, ridgeline KDE)
- `run_settings.json` capturing all parameters for reproducibility

### Export History

Browse past export runs, view their settings, and:

- Reload settings from a previous run into the Setup tab.
- Delete old export directories.
- Download any export run as a ZIP file.

### Plot Gallery

Browse exported plots organized by peptide and reference. Download individual
plots without needing to navigate the filesystem.

See [export/layout.md](export/layout.md) for the full directory structure
and [export/formats.md](export/formats.md) for file schemas.

---

## Tips

- **Peptide preview** -- Use the "Preview peptides in idXML" button in Setup
  before searching. This shows a sample of peptides in your files and helps
  catch typos.
- **Same-name files** -- If you select files with identical names from different
  directories, MSpecViz automatically disambiguates them using parent directory
  components.
- **Re-runs with same settings** -- MSpecViz detects matching export directories
  by settings fingerprint. Existing plots are skipped, saving time on re-exports.
- **Large datasets** -- Increase the worker count in Performance Settings to
  speed up export. Each worker runs in a separate process to avoid GIL
  contention.
- **Koina fallback** -- Enable Koina AI prediction when you lack experimental
  reference spectra. The predicted spectrum serves as a reasonable baseline
  for mirror plot comparison, though experimental references are always
  preferred.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| No files shown in the directory browser | Verify the path exists and you have read permissions. Check the `data_directories.default` setting in `streamlit_app/config/settings.json`. |
| "No matches found" after peptide search | Verify the peptide sequences exist in your idXML files using the Preview button. Sequences are case-insensitive but must match the base (unmodified) sequence. |
| PSMs not resolved to mzML files | Ensure the mzML files referenced in the idXML metadata are among your selected files. The filenames must match. |
| openms_insight components not rendering | Confirm openms_insight 0.1.11 is installed. Restart the Streamlit server if components fail to load. |
| Koina prediction returns no results | Check network connectivity (Koina uses a remote API). Some models only support specific peptide lengths or charge states. |
| Slow first run | Normal for large mzML files. Subsequent runs with the same settings reuse existing export directories. |
| Export fails with permission error | Verify write permissions to the export directory. Check available disk space. |
| Worker processes crash silently | Re-run with `workers: 1` in the YAML (CLI) or set the workers slider to 1 in the Setup tab to surface the underlying traceback. |
| Stale results after fixing a bug | Delete `<workspace>/.mspecviz_cache/`. The app does not invalidate its own Parquet cache automatically. See [settings/app.md](settings/app.md) for the cache key in `settings.json`. |
