# CLI Pipeline Reference

MSpecViz includes a terminal-based pipeline (`cli.py`) that runs the full analysis without opening a browser. It reuses the same backend as the Streamlit app, so results are identical.

## Quick Start

```bash
# 1. Ensure you have activated the correct conda environment
conda activate mspecviz

# 2a. Minimal run — just a config file
python cli.py --config examples/config_cli.yaml

# 2b. With verbose logging
python cli.py --config examples/config_cli.yaml --verbose

# 2c. Override tolerance and add neutral losses
python cli.py --config examples/config_cli.yaml --mz-tolerance 0.05 --add-losses
```

An example configuration file with all available options is provided at [`examples/config_cli.yaml`](../examples/config_cli.yaml).

## Configuration Priority

Settings are merged from three layers (highest priority wins):

| Priority | Source | Example |
|----------|--------|---------|
| 1 (highest) | CLI arguments | `--mz-tolerance 0.05` |
| 2 | YAML config file | `annotation: mz_tolerance: 0.02` |
| 3 (lowest) | `settings.json` | App-wide defaults |

This means you can set reasonable defaults in the YAML and override individual values on the command line for one-off experiments.

## Configuration File

The YAML config has seven sections:

### `input` — File paths

```yaml
input:
  idxml_files:
    - /path/to/data/*.idXML          # glob patterns supported
    - /path/to/more/**/*.idXML       # recursive with **
  mzml_files:
    - /path/to/data/*.mzML
  # parquet_file: /path/to/cache.parquet  # skip search + mzML loading
```

- **Glob patterns** are expanded before processing (`*`, `**`, `?` all work).
- **`parquet_file`** (optional): If set, the pipeline skips peptide search and mzML loading. Use this for instant re-runs after the first analysis.

### `peptides` — Target sequences

```yaml
# Option 1: File path (one peptide per line, recommended)
peptides: examples/peptides.txt

# Option 2: Extract all peptides from idXML files
peptides: all

# Option 3: Inline list (for quick tests)
peptides: [LGGNEQVTR, GAGSSEPVTGLDAK]
```

Modifications like `(Oxidation)` are automatically stripped to base sequences. Duplicates are removed.

### `reference` — Reference file identification

```yaml
reference:
  substring: "Synpep"    # case-insensitive match against mzML filenames
```

mzML files whose name contains the substring are treated as reference spectra. All others become queries. Mirror plots pair each query with its matching reference(s).

Omit or leave empty for single-spectrum mode (no mirroring).

### `annotation` — Peak annotation

```yaml
annotation:
  ion_types: [b, y]              # a, b, c, x, y, z, precursor, immonium
  mz_tolerance: 0.02             # Da
  add_losses: false               # H2O and NH3 neutral losses
  add_first_prefix_ion: true
  remove_precursor_peaks: false
  intensity_mode: Relative        # "Relative" (0-1) or "Raw"
```

### `export` — Output settings

```yaml
export:
  directory: ./exports
  project_name: my_analysis
  formats: [png]                  # [png], [svg], or [png, svg]
  dpi: 300
  figsize: [12, 8]                # width, height in inches
  plot_unpaired: false            # plot spectra without matching reference
```

### `performance` — Parallelism

```yaml
performance:
  workers: 4   # number of parallel processes for plot generation (1 = sequential)
```

### `cache` — Parquet auto-caching

```yaml
cache:
  enabled: true
  subdirectory: .mspecviz_cache   # relative to export directory
```

After the first run the loaded spectra are cached as Parquet. Subsequent runs with identical inputs load from cache automatically. See [Caching System](caching.md) for details.

## CLI Arguments

Every config setting can be overridden from the command line.

### Input overrides

| Flag | Description |
|------|-------------|
| `--config PATH` | **(required)** Path to YAML configuration file |
| `--idxml FILE [FILE ...]` | idXML file paths (overrides config) |
| `--mzml FILE [FILE ...]` | mzML file paths (overrides config) |
| `--parquet PATH` | Parquet file (skips search + mzML loading) |
| `--peptides PATH` | Peptides .txt file (one per line) |
| `--reference-substring STR` | Reference file identifier substring |

### Annotation overrides

| Flag | Description |
|------|-------------|
| `--ion-types TYPE [TYPE ...]` | Ion types (e.g., `b y precursor`) |
| `--mz-tolerance FLOAT` | m/z tolerance in Da |
| `--add-losses` | Enable neutral losses (H2O, NH3) |
| `--no-first-prefix-ion` | Disable first prefix ion |
| `--remove-precursor-peaks` | Remove precursor peaks before annotation |
| `--raw-intensity` | Use raw intensities instead of normalized |

### Export overrides

| Flag | Description |
|------|-------------|
| `--export-dir PATH` | Export directory |
| `--project-name NAME` | Project name for export folder |
| `--formats FMT [FMT ...]` | Export formats (`png`, `svg`) |
| `--dpi INT` | DPI for image export |
| `--plot-unpaired` | Plot spectra without matching reference |

### Other flags

| Flag | Description |
|------|-------------|
| `--no-cache` | Disable Parquet auto-caching |
| `--workers INT` | Number of parallel worker processes for plot generation (default: 4) |
| `-v`, `--verbose` | Verbose output with debug details |
| `-q`, `--quiet` | Suppress all output except errors |

## Pipeline Phases

The CLI runs seven sequential phases. In verbose mode (`--verbose`), each phase prints detailed information.

### Phase 0: Input validation

Checks that all specified idXML/mzML files exist. If `peptides: all`, extracts every unique peptide from the idXML files.

### Phase 1: Auto-cache detection

When caching is enabled, computes a fingerprint from the input files + peptides + reference substring. If a matching Parquet cache exists, switches to the fast Parquet path automatically.

### Phase 2: Peptide search

Searches the idXML files for PSMs matching the target peptides. Reports total PSMs found and the number of peptide variants.

### Phase 3: Spectrum loading

Streams spectra from mzML files (or loads instantly from Parquet). Reference spectra are loaded first, then queries are processed one-by-one as they arrive.

### Phase 4a: Annotation

For each query spectrum arriving from the stream:
1. Generates theoretical fragment ions for the peptide sequence
2. Matches experimental peaks against theoretical within the m/z tolerance
3. Collects the annotated spectrum data into a task list for plotting

### Phase 4b: Parallel plot generation

The collected tasks are processed in parallel using multiple worker processes (`--workers`, default 4):
1. Each worker receives one annotated spectrum pair and independently computes metrics, creates a mirror plot (or single-spectrum plot), and saves it to disk
2. Workers run in separate processes with their own matplotlib backends, so memory is released after each plot
3. A progress bar tracks completion as tasks finish
4. If `--workers 1`, processing falls back to sequential mode (useful for debugging)

### Phase 5: Summary export

Writes a comparison summary as CSV and XLSX with per-pair metrics: R² (matched/binned), cosine similarity, spectral entropy, Dice, Jaccard, and shared ion ratio.

### Phase 6: Cache save

If caching is enabled and this was a fresh run (not loaded from cache), saves the loaded spectra as a Parquet file for instant re-use.

## Output Structure

```
exports/
└── my_analysis_20260220_143000/
    ├── png/
    │   ├── LGGNEQVTR/
    │   │   ├── LGGNEQVTR_sample1_vs_synpep1.png
    │   │   └── LGGNEQVTR_sample2_vs_synpep1.png
    │   └── GAGSSEPVTGLDAK/
    │       └── ...
    ├── svg/                          # if formats includes svg
    │   └── ...
    ├── comparison_summary.csv
    ├── comparison_summary.xlsx
    └── .mspecviz_cache/              # auto-cache (if enabled)
        └── <fingerprint>.parquet
```

## Examples

### Basic analysis

```bash
python cli.py --config examples/config_cli.yaml
```

### Verbose output for debugging

```bash
python cli.py --config examples/config_cli.yaml --verbose
```

### Override peptides and tolerance

```bash
python cli.py --config examples/config_cli.yaml \
  --peptides my_peptides.txt \
  --mz-tolerance 0.05
```

### Export both PNG and SVG

```bash
python cli.py --config examples/config_cli.yaml --formats png svg
```

### All peptides, no cache

```bash
python cli.py --config examples/config_cli.yaml \
  --peptides all \
  --no-cache
```

### Re-run from Parquet cache

```bash
python cli.py --config examples/config_cli.yaml \
  --parquet exports/my_analysis_20260220_143000/.mspecviz_cache/abc123.parquet
```

### Quiet mode for scripting

```bash
python cli.py --config examples/config_cli.yaml --quiet
# Exit code: 0 = success, 1 = processing failure, 2 = config error
```

---

## Standalone Scripts

MSpecViz includes two additional standalone scripts beyond the main CLI pipeline:

### Summary & Per-Peptide Plots (`src/summary_plots.py`)

Generates cross-peptidoform summary plots and per-peptide dumbbell plots from a `comparison_summary.csv` produced by the main pipeline.

```bash
# Generate all summary plots (dual-axis, ridgeline, heatmap)
python -m src.summary_plots exports/my_run/comparison_summary.csv

# Specify output directory and metric
python -m src.summary_plots exports/my_run/comparison_summary.csv \
  -o exports/my_run -m "R2 Binned"

# Generate only specific plot types
python -m src.summary_plots exports/my_run/comparison_summary.csv \
  --plots dual_axis heatmap

# Generate per-peptide dumbbell plots (R² binned vs R² match per query)
python -m src.summary_plots exports/my_run/comparison_summary.csv \
  --peptide-detail

# Generate dumbbell plots for specific peptidoforms only
python -m src.summary_plots exports/my_run/comparison_summary.csv \
  --peptide-detail --peptidoforms "LGGNEQVTR/2" "GAGSSEPVTGLDAK/3"

# Export as SVG
python -m src.summary_plots exports/my_run/comparison_summary.csv \
  --formats svg --dpi 300

# List available metrics in the CSV
python -m src.summary_plots exports/my_run/comparison_summary.csv --list-metrics
```

| Flag | Description |
|------|-------------|
| `csv` | **(positional, required)** Path to `comparison_summary.csv` |
| `-o`, `--output-dir` | Output directory (default: same as CSV) |
| `-m`, `--metric` | Metric to visualize (default: `R2 Binned`) |
| `--plots` | Summary plot types: `dual_axis`, `ridgeline`, `heatmap` (default: all) |
| `--formats` | Export formats: `png`, `svg` (default: `png`) |
| `--dpi` | Image DPI (default: 300) |
| `--title` | Title prefix for plots |
| `--list-metrics` | List available metrics in the CSV and exit |
| `--peptide-detail` | Generate per-peptide dumbbell plots instead of summary plots |
| `--peptidoforms` | Specific peptidoforms to plot (default: all) |

**Summary plots** are saved to `{output_dir}/summary/`. **Dumbbell plots** are saved alongside mirror plots in `{output_dir}/png/{Peptide}/` (and/or `svg/`).

### Parquet Conversion Tool (`tools/convert_mzml_to_parquet.py`)

Converts idXML + mzML files into a single Parquet file for instant loading on subsequent runs. Useful for large datasets or batch processing.

```bash
python tools/convert_mzml_to_parquet.py \
  --idxml /path/to/*.idXML \
  --mzml-dir /path/to/mzml_files/ \
  --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
  --output exports/my_spectra.parquet \
  --reference-substring "Synpep"
```

| Flag | Description |
|------|-------------|
| `--idxml` | **(required)** Path(s) to idXML file(s) |
| `--mzml-dir` | **(required)** Directory containing mzML files |
| `--peptides` | **(required)** Comma-separated list of target peptide sequences |
| `--output` | **(required)** Output Parquet file path |
| `--reference-substring` | Substring to identify reference mzML files |
| `-v`, `--verbose` | Enable verbose logging |

See [Parquet Conversion Tool](../tools/README.md) for full documentation.
