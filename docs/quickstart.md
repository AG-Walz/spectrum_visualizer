# Quick Start

An end-to-end run with the bundled test data, plus the CLI equivalent. Deeper
references: [architecture.md](architecture.md), [settings.md](settings.md),
[cli.md](cli.md), [metrics.md](metrics.md), [export.md](export.md).

---

## What you need

Two file types from the **same experiment**: `.idXML` (peptide-spectrum matches from
a search engine) and `.mzML` (raw MS/MS spectra). The idXML references spectra inside
the mzML by spectrum ID.

Test data for a verification run (expected: 2 PSMs — ASDDFDQASY, RTEFNLNQY — 373 peaks):

| Input | Path |
|---|---|
| idXML | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/ProT001_1_pout_filtered.idXML` |
| mzML dir | `/mnt/volume/nacho/data/ProCa_runs/small_dataset/` |
| Peptide list | `/mnt/volume/nacho/data/peptide_lists/peptides.txt` |

---

## Run it

**Streamlit** (run from `streamlit_app/` so `.streamlit/config.toml` is picked up;
workspaces still resolve to the project root, so `cd` is safe):

```bash
conda activate mspecviz
cd streamlit_app && streamlit run app.py     # → http://localhost:8501
```

**CLI** (headless, from the project root):

```bash
python -m cli cli/config/config_example.yaml
```

YAML-only — copy the example and edit it for per-run tweaks; no flag overrides.
`MSPECVIZ_LOG_LEVEL=ERROR` for a quiet run. Full reference: [cli.md](cli.md).

---

## The four tabs

1. **Setup** — select idXML + mzML files (directory browsers), enter target base
   sequences (case-insensitive; modified forms are found automatically), set the
   reference substring, and the annotation / Koina / performance options. Defaults
   live in `core/defaults.py`; see [settings.md](settings.md).
2. **Load & Annotate** — *Look Up Peptides* groups PSMs by peptidoform and resolves
   their mzML; *Load & Annotate Spectra* streams the spectra, annotates fragments,
   pairs queries with references (or a Koina prediction), computes the 7 metrics, and
   builds the openms_insight DataFrames.
3. **Explore** — three sub-views: the linked openms_insight Interactive Explorer
   (Table + SequenceView + LinePlot + Heatmap), Plotly Mirror Plots (query up,
   reference down), and the Comparison Metrics tables. See [metrics.md](metrics.md).
4. **Export** — versioned output directory with Parquet, CSV/XLSX, plots (PNG/SVG/PDF/
   HTML), summary plots, and `run_settings.json`, plus a plot gallery and run history.
   See [export.md](export.md).

**Reference substring.** A filename substring marking reference mzMLs (e.g. `Synpep`)
enables mirror plots + metrics. Leave it empty for single-spectrum mode (annotate
only, no comparison). When no experimental reference exists and Koina is enabled, an
AI-predicted reference is generated as a fallback.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| No files in the browser | Check the path exists and `data_directories.default` in `settings.json`. |
| "No matches" after lookup | Verify the sequences exist via the Preview button; match the base (unmodified) sequence. |
| PSMs not resolved to mzML | The mzML files referenced in the idXML must be among your selection (filenames must match). |
| openms_insight components blank | Confirm `openms_insight 0.1.11` is installed; restart the server. |
| Koina returns nothing | Check connectivity; some models only support specific lengths/charges. |
| Workers crash silently | Re-run with `workers: 1` (YAML) or the slider at 1 to surface the traceback. |
| Stale results after a bug fix | Delete `<workspace>/.mspecviz_cache/` — the app never invalidates its own Parquet cache. |
