# MSpecViz

MSpecViz is a peptide-spectrum match (PSM) analysis and visualization toolkit.
The core pipeline is pure Python and ships with two interchangeable front ends:

- **Streamlit UI** — interactive 4-tab workflow built on the
  [openms_insight](https://pypi.org/project/openms-insight/) library with
  Plotly mirror plots and quantitative comparison metrics.
- **CLI** — headless pipeline driven by a YAML config, suitable for batch jobs
  and servers without a browser.

Both front ends share the same `core/` core (peptide lookup, annotation,
metrics, export). Pick the interface that fits your workflow; neither is required.

---

## Table of Contents

### Reading guides (in `docs/`)

Start with the **Quick Start**; the rest are topic references. They describe
structure and rationale, pointing at the code (`core/defaults.py`, the schema
constants) for anything that would otherwise drift.

- [Quick Start](docs/quickstart.md) — end-to-end run with test data + the CLI equivalent
- [Architecture](docs/architecture.md) — three layers, pipeline flow, house rules, the `core/` map, front ends, tests
- [Settings](docs/settings.md) — sources, precedence, `settings.json`, workspaces, session-state contract
- [CLI](docs/cli.md) — invocation, YAML config, recipes, debugging
- [Export](docs/export.md) — versioned dirs, Parquet/CSV formats, plots, gallery & history
- [Metrics](docs/metrics.md) — the 7 comparison metrics with formulas + how to choose

### Jump to a README section

- [Key Features](#key-features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Configuration at a glance](#configuration-at-a-glance)
- [Running Tests](#running-tests)
- [Troubleshooting](#troubleshooting)

---

## Key Features

- **Interactive Explorer** — openms_insight Table, LinePlot, SequenceView, and
  Heatmap components linked through a shared state manager. Click a PSM row to
  instantly view its annotated spectrum and sequence coverage.
- **Mirror Plots** — Plotly interactive and Matplotlib static mirror plots
  comparing query spectra against experimental or AI-predicted references.
- **7 Comparison Metrics** — R² matching, R² binned, cosine similarity,
  spectral entropy similarity, Dice coefficient, Jaccard index, and shared
  peak ratio. See [docs/metrics.md](docs/metrics.md).
- **Koina AI Prediction** — 14 deep-learning models (MS2PIP, Prosit,
  AlphaPept, UniSpec, PredFull) for generating predicted reference spectra
  when no experimental reference is available.
- **Parallel Export** — `ProcessPoolExecutor`-based batch generation of PNG,
  SVG, PDF, and HTML plots with versioned export directories and ZIP
  download.
- **Summary Plots** — Dual-axis bar/box and ridgeline KDE visualizations of
  metric distributions across peptidoforms.

---

## Prerequisites

- Python 3.11
- Conda (Miniconda or Anaconda)
- The `mspecviz` conda environment

Key dependencies (installed in the environment):

`streamlit >= 1.48`, `openms_insight 0.1.11`, `pyopenms`, `polars`, `plotly`,
`matplotlib`, `ms_entropy`, `koinapy`, `psm_utils`, `scipy`, `PyYAML`.

---

## Installation

Clone or copy the project, then activate the conda environment:

```bash
conda activate mspecviz
```

The environment is located at `/mnt/volume/bin/miniconda/envs/mspecviz/`.

To recreate it from scratch:

```bash
conda env create -f environment.yml
```

---

## Quick Start

MSpecViz has two interfaces — pick whichever fits your workflow.

### Streamlit (interactive UI)

Run from the `streamlit_app/` directory so that `.streamlit/config.toml` is
picked up:

```bash
cd streamlit_app
streamlit run app.py
```

Or with the full path to the environment's Streamlit:

```bash
cd streamlit_app
/mnt/volume/bin/miniconda/envs/mspecviz/bin/streamlit run app.py
```

Open http://localhost:8501 in your browser. For a remote server:

```bash
ssh -L 8501:localhost:8501 user@server
```

Workspaces, exports, and workspace defaults all resolve relative to the
**project root**, so you can safely `cd` into `streamlit_app/` without
reshuffling data directories.

### CLI (headless pipeline)

1. Copy the bundled example to a run-specific file and edit the paths at
   the top (inputs, peptides, reference substring, export directory):

   ```bash
   cp cli/config/config_example.yaml cli/config/my_run.yaml
   ```

   The example lists **every** field the CLI understands, grouped by how
   often you need to touch it (required inputs first, sensible defaults in
   the middle, Koina / advanced knobs at the bottom). Most runs only need
   the first block.

2. Run from the project root:

   ```bash
   python -m cli cli/config/my_run.yaml
   ```

The CLI is YAML-only — there are no command-line flag overrides. For
per-run tweaks (different workers, tighter `mz_tolerance`, etc.), copy
the YAML and edit it. Quiet runs: `MSPECVIZ_LOG_LEVEL=ERROR python -m
cli …`. See [docs/cli.md](docs/cli.md) for the full CLI reference, or
[docs/quickstart.md](docs/quickstart.md) for a walkthrough with the bundled
test data.

---

## Project Structure

Three top-level packages keep core logic separable from the UI runners:
`core/` (pure Python — pipeline, helpers, every tunable default), and two
sibling front ends that both import it, `streamlit_app/` (the 4-tab UI) and
`cli/` (the headless YAML-driven runner). `tests/` mirrors `core/` and covers
all three.

**Rules the codebase follows:**

- `core/` must never import from `streamlit_app/` or `cli/` — the core stays
  testable and reusable.
- Both front ends depend on `core/`; anything shared goes there (usually
  `core/pipeline/` for orchestration or `core/helpers.py` for pure helpers).
- No `.py` file exceeds 250 lines; no folder holds more than 10 scripts —
  past either cap, split into a subpackage with a re-exporting `__init__.py`.

The layered design, pipeline flow, and a tour of the `core/` subpackages and
front ends live in [docs/architecture.md](docs/architecture.md).

---

## Configuration at a glance

MSpecViz reads configuration from a small set of sources. Full reference
with every key documented is in
[docs/settings.md](docs/settings.md); here is the map:

| File | Scope | What it controls |
|---|---|---|
| `core/defaults.py` (`DEFAULTS`) | Single source of truth for tunable defaults | Every parameter's fallback value (ion types, tolerances, export DPI, Koina settings, …) |
| `streamlit_app/config/settings.json` | App-wide environmental config | Koina model catalogue, allowed ion types, data directory, performance tuning, ion color palette |
| `workspaces/<name>/params.json` | One per workspace | Sparse record of user-edited widget values; missing keys fall back to `DEFAULTS` |
| `cli/config/*.yaml` | One per CLI run | YAML config combining input files, peptides, annotation, export, and Koina sections |
| `MSPECVIZ_LOG_LEVEL` env var | Single invocation | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) — the only out-of-band override |

**Precedence for CLI runs:** `DEFAULTS` < YAML. The CLI takes one positional
argument (the YAML path); there are no flag overrides — edit the YAML for
per-run tweaks.

**Precedence for Streamlit runs:** `DEFAULTS` < workspace `params.json`
< live widget state in session. (`settings.json` is read in parallel for
environmental config and enumerations; it does not supply tunable
defaults.)

---

## Running Tests

```bash
conda run -n mspecviz python -m pytest tests/ -v
```

Test groups, mirroring the source layout:

- `tests/core/` — unit + integration tests for the core (grouped into
  subfolders mirroring the `core/` subpackages)
- `tests/cli/` — CLI integration tests (invoke `python -m cli` as a
  subprocess)
- `tests/streamlit/` — GUI smoke tests (via `streamlit.testing.v1.AppTest`)
  plus export-pipeline unit tests

Marker-driven selection:

```bash
pytest tests/ -v -m "not integration"    # skip the slower end-to-end tests
pytest tests/core/ -k "annotation or metrics" -v   # annotation + metrics units
```

---

## Troubleshooting

See the "Troubleshooting" section of
[docs/quickstart.md](docs/quickstart.md#troubleshooting) for the common
failure modes (missing files, empty peptide matches, Koina connectivity,
permission errors on export).

For deeper issues:

- Streamlit components render blank → confirm `openms_insight 0.1.11` is
  installed (`pip show openms_insight`) and restart the server.
- Stale cache after fixing a loading bug → delete the workspace's
  `.mspecviz_cache/` directory. The app never invalidates its own Parquet
  cache automatically.
- Worker processes crash silently → set `workers: 1` in your CLI YAML or
  drag the workers slider to 1 in the Setup tab to surface the real error.
