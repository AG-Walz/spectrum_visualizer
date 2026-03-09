# MSpecViz Quick Start Guide

Get up and running in 5 minutes. For detailed explanations, see the [Full Documentation](DOCUMENTATION.md).

---

## Installation

```bash
conda env create -f environment.yml   # first time only
conda activate mspecviz
```

**Streamlit app** (interactive, browser-based):
```bash
streamlit run app.py
```
Open http://localhost:8501. If connecting to a remote server, see [Remote Access](../README.md#remote-access-ssh-tunnel).

**CLI pipeline** (batch processing, no browser needed):
```bash
python cli.py --config examples/config_cli.yaml
```
See the [CLI Reference](cli.md) for the full configuration reference.

---

## What You Need

| File Type | Extension | Contains |
|-----------|-----------|----------|
| **Identification files** | `.idXML` | Peptide-spectrum matches (PSMs) from a search engine |
| **Spectrum files** | `.mzML` | Raw MS/MS spectra (m/z + intensity arrays) |

Both file types must come from the same experiment. The idXML files contain references to spectra inside the mzML files.

---

## The 7-Step Workflow

MSpecViz is organized into 7 sequential tabs. Work through them left to right.

### Step 1: Files

Select your input files using the built-in directory browser.

- **Left column**: Select one or more `.idXML` files
- **Right column**: Select one or more `.mzML` files
- Navigate directories by clicking folder buttons or typing a path
- Select files from multiple directories — selections accumulate
- Use the filter box to narrow down long file lists

**Parquet shortcut**: If you have a previously converted `.parquet` file, paste its path at the bottom to skip directly to Step 6 (Process).

### Step 2: Peptides

Enter target peptide sequences, one per line:

```
LGGNEQVTR
GAGSSEPVTGLDAK
VEATFGVDESNAK
```

- Enter **base sequences** (unmodified) — modified variants are found automatically
- Case insensitive
- Click "Preview in idXML" to see a sample of peptides present in your files

### Step 3: Reference

Enter a substring that identifies your reference mzML files (e.g., `Synpep`, `Standard`, `Ref`).

- Files whose names contain the substring are categorized as **Reference**
- All other files are categorized as **Query**
- A preview shows how your files will be split

**If you leave this blank**, every spectrum is plotted individually with annotations (no mirror plots).

### Step 4: Settings

Configure analysis parameters:

- **Ion types**: Which fragment ions to annotate (default: b, y)
- **m/z tolerance**: How close an experimental peak must be to a theoretical fragment to count as a match (default: 0.02 Da)
- **Neutral losses**: Include H₂O and NH₃ loss peaks
- **First prefix ion**: Include a₁/b₁ ions
- **Remove precursor peaks**: Strip precursor ions from plots
- **Intensity mode**: Relative (normalized to base peak) or Raw

### Step 5: Search

Click **Search Peptides**. MSpecViz searches your idXML files for the target peptides, finds all PSMs, and resolves which mzML file contains each spectrum.

The results table shows each peptide variant, its modification state, how many PSMs were found, and whether the mzML files were successfully resolved.

### Step 6: Process

Configure export settings (project name, directory, formats) and click **Process & Export**.

MSpecViz will:
1. Load reference spectra from mzML (or Parquet cache)
2. Load each query spectrum
3. Annotate peaks with theoretical fragment ions
4. Generate mirror plots (query on top, reference on bottom)
5. Compute comparison metrics (cosine similarity, R², entropy, etc.)
6. Export plots as PNG (and optionally SVG)
7. Write a comparison summary spreadsheet (CSV + XLSX)

Progress bars show reference loading and query processing separately, with time estimates.

### Step 7: Results

View a summary of what was exported: number of plots, file counts, export location, and any errors. Expandable sections list individual exported files.

---

## Quick Tips

- **Same-name files**: If you select files with identical names from different directories, MSpecViz automatically disambiguates them by adding parent directory components.
- **Auto-caching**: After the first successful run, MSpecViz saves a Parquet cache. Subsequent runs with the same inputs load instantly. See [Caching](caching.md).
- **Unpaired spectra**: By default, spectra without a matching reference are skipped. Check "Plot unpaired spectra" in Export Settings to include them as single-spectrum plots.
- **Large datasets**: For very large datasets, use the [CLI conversion tool](../tools/README.md) to pre-convert to Parquet once, then load instantly every time.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No files shown in browser | Check the path exists and you have read permissions |
| "No matches found" after search | Verify the peptide sequences exist in your idXML files (use Preview) |
| PSMs could not be mapped to mzML | Ensure the mzML files listed in the idXML metadata are among your selected files |
| Slow first run | Normal for large mzML files. Subsequent runs use the Parquet cache |
| White screen flash on button click | This should not happen. If it does, restart the Streamlit server |

---

## Further Reading

| Topic | Link |
|-------|------|
| Full user guide | [DOCUMENTATION.md](DOCUMENTATION.md) |
| Caching system | [caching.md](caching.md) |
| Comparison metrics | [spectrum_metrics.md](spectrum_metrics.md) |
| File selection details | [file_selection.md](file_selection.md) |
| Annotation settings | [annotation_settings.md](annotation_settings.md) |
| Mirror plot interpretation | [mirror_plots.md](mirror_plots.md) |
| Export & output formats | [export.md](export.md) |
| Parquet conversion tool | [tools/README.md](../tools/README.md) |

---

## References

MSpecViz is built on the [OpenMS WebApps template](https://abi-services.cs.uni-tuebingen.de/streamlit-template/):

Muller, T. D., Siraj, A., et al. (2025). OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. *Journal of Proteome Research*. [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872)
