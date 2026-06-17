# MSpecViz

MSpecViz is a peptide-spectrum match (PSM) analysis and visualization toolkit.
The core pipeline is pure Python and ships with two interchangeable front ends:

- **Streamlit UI** — interactive 4-tab workflow built on the
  [openms_insight](https://pypi.org/project/openms-insight/) library with
  Plotly mirror plots and quantitative comparison metrics.
- **CLI** — headless pipeline driven by a YAML config, suitable for batch jobs
  and servers without a browser.

Both front ends share the same `core/` package (peptide lookup, annotation,
metrics, export). Pick the interface that fits your workflow; neither is required.

## Setup

Requires Python 3.11 and Conda. Create and activate the environment (the create
step is first-run only):

```bash
conda env create -f environment.yml
conda activate mspecviz
```

## Run it

**Streamlit** (run from `streamlit_app/` so `.streamlit/config.toml` is picked up;
workspaces and exports still resolve to the project root, so the `cd` is safe):

```bash
cd streamlit_app && streamlit run app.py     # → http://localhost:8501
```

**CLI** (headless, from the project root) — copy the bundled example, edit its paths,
and run. It is YAML-only with no flag overrides; `MSPECVIZ_LOG_LEVEL=ERROR` gives a
quiet run. Full reference: [CLI → YAML config](docs/cli.md#yaml-config).

```bash
cp cli/config/config_example.yaml cli/config/my_run.yaml
python -m cli cli/config/my_run.yaml
```

You need two file types from the **same experiment**: `.idXML` (peptide-spectrum
matches from a search engine) and `.mzML` (raw MS/MS spectra); the idXML references
spectra inside the mzML by spectrum ID. Bundled test data for a verification run
(expected: 2 PSMs — ASDDFDQASY, RTEFNLNQY — 373 peaks):

| Input | Path |
|---|---|
| idXML | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/ProT001_1_pout_filtered.idXML` |
| mzML dir | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/` |
| Peptide list | `/mnt/volume/nacho/data/peptide_lists/peptides.txt` |

## The four tabs (Streamlit)

1. **Setup** — select idXML + mzML files (directory browsers), enter target base
   sequences (case-insensitive; modified forms are found automatically), set the
   reference substring, and the annotation / Koina / performance options. Defaults
   live in `core/defaults.py`; see [Settings → sources and precedence](docs/settings.md#sources-and-precedence).
2. **Load & Annotate** — *Look Up Peptides* groups PSMs by peptidoform and resolves
   their mzML; *Load & Annotate Spectra* streams the spectra, annotates fragments,
   pairs queries with references (or a Koina prediction), computes the 7 metrics, and
   builds the openms_insight DataFrames.
3. **Explore** — three sub-views: the linked openms_insight Interactive Explorer
   (Table + SequenceView + LinePlot + Heatmap), Plotly Mirror Plots (query up,
   reference down), and the Comparison Metrics tables. See
   [Metrics → the seven metrics](docs/metrics.md#the-seven-metrics).
4. **Export** — versioned output directory with Parquet, CSV/XLSX, plots (PNG/SVG/PDF/
   HTML), summary plots, and `run_settings.json`, plus a plot gallery and run history.
   See [Export → directory layout and versioning](docs/export.md#directory-layout-and-versioning).

**Reference substring.** A filename substring marking reference mzMLs (e.g. `Synpep`)
enables mirror plots + metrics. Leave it empty for single-spectrum mode (annotate
only, no comparison). When no experimental reference exists and Koina is enabled, an
AI-predicted reference is generated as a fallback. Full matching rules:
[CLI → reference substring and single-spectrum mode](docs/cli.md#reference-substring-and-single-spectrum-mode).

## Documentation

Topic references — structure and rationale, pointing at the code
(`core/defaults.py`, the schema constants) for anything that would otherwise drift.

- [Architecture](docs/architecture.md) — three layers, pipeline flow, house rules, the `core/` map
- [Settings](docs/settings.md) — sources, precedence, `settings.json`, workspaces, session-state contract
- [CLI](docs/cli.md) — invocation, YAML config, recipes, debugging
- [Export](docs/export.md) — versioned dirs, Parquet/CSV formats, plots, gallery & history
- [Metrics](docs/metrics.md) — the 7 comparison metrics with formulas + how to choose

## Troubleshooting

| Problem | Fix |
|---|---|
| No files in the browser | Check the path exists and `data_directories.default` in `settings.json`. |
| "No matches" after lookup | Verify the sequences exist via the Preview button; match the base (unmodified) sequence. |
| PSMs not resolved to mzML | The mzML files referenced in the idXML must be among your selection (filenames must match). |
| openms_insight components blank | Confirm `openms_insight 0.1.13` is installed; restart the server. |
| Koina returns nothing | Check connectivity; some models only support specific lengths/charges. |
| Workers crash silently | Re-run with `workers: 1` (YAML) or the slider at 1 to surface the traceback. |
| Stale results after a bug fix | Delete `<workspace>/.mspecviz_cache/` — the app never invalidates its own Parquet cache. |

## Running tests

```bash
conda run -n mspecviz python -m pytest tests/ -v
pytest tests/ -m "not integration"    # skip the slower end-to-end tests
```
