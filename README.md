# MSpecViz

MSpecViz is a peptide-spectrum match (PSM) analysis and visualization toolkit.
The core pipeline is pure Python and ships with two interchangeable front ends:

- **Streamlit UI** — interactive 4-tab workflow built on the
  [openms_insight](https://pypi.org/project/openms-insight/) library with
  Plotly mirror plots and quantitative comparison metrics.
- **CLI** — headless pipeline driven by a YAML config, suitable for batch jobs
  and servers without a browser.

Both front ends share the same `mspecviz/` core (search, annotation, metrics,
export). Pick the interface that fits your workflow; neither is required.

---

## Table of Contents

### Reading guides (in `docs/`)

Start with the **Quick Start** for an end-to-end walkthrough; the other
guides are reference material organised by topic, with each file kept under
~200 lines for easy skimming.

**Getting started**
- [Quick Start](docs/quickstart.md) — 4-tab Streamlit walkthrough + the `python -m cli` equivalent

**Architecture**
- [Overview](docs/architecture/overview.md) — 3-layer layout, pipeline flow, house rules
- [Core package](docs/architecture/core.md) — tour of every `mspecviz/` subpackage
- [Front-ends](docs/architecture/frontends.md) — `streamlit_app/`, `cli/`, `tests/`

**Settings**
- [Overview](docs/settings/overview.md) — precedence rules + quick-lookup table
- [`settings.json`](docs/settings/app.md) — every app-wide key explained
- [Workspace + session state](docs/settings/workspace.md) — `params.json`, `peptide_mirror.*` keys, `DEFAULTS` fallbacks
- [CLI YAML + flags](docs/settings/cli.md) — every CLI field and command-line flag

**CLI**
- [Usage](docs/cli/usage.md) — invocation, minimum config, exit codes, output
- [Recipes & debugging](docs/cli/recipes.md) — peptide spec, common recipes, troubleshooting

**Metrics**
- [7 comparison metrics](docs/metrics/comparison.md) — R2m, R2b, Cos, Ent, Dice, Jacc, ShR with formulas
- [17 ion match statistics](docs/metrics/ion-match-stats.md) — `ims_*` columns
- [Choosing and interpreting](docs/metrics/interpretation.md) — decision guide + common patterns

**Export**
- [Directory layout](docs/export/layout.md) — versioned dirs + `run_settings.json`
- [Data file formats](docs/export/formats.md) — Parquet schemas + CSV/XLSX columns
- [Plots, gallery, history](docs/export/plots.md) — mirror + summary plots, Streamlit browsers

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
  peak ratio. See [docs/metrics/comparison.md](docs/metrics/comparison.md).
- **17 Ion Match Statistics** — PyOpenMS `SpectrumAnnotator` metrics including
  coverage, annotated-ion fraction, fragment mass errors, and series analysis.
- **Koina AI Prediction** — 14+ deep-learning models (MS2PIP, Prosit,
  AlphaPept, UniSpec, PredFull) for generating predicted reference spectra
  when no experimental reference is available.
- **Parallel Export** — `ProcessPoolExecutor`-based batch generation of PNG,
  SVG, PDF, and HTML plots with versioned export directories and ZIP
  download.
- **Summary Plots** — Dual-axis bar/box, ridgeline KDE, and multi-metric
  heatmap visualizations of metric distributions across peptidoforms.

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
   python -m cli -c cli/config/my_run.yaml
   ```

Any field can be overridden on a single invocation with a flag — for
example `--workers 8 --mz-tolerance 0.05`. See the comment block at the
bottom of the example YAML, or [docs/cli/usage.md](docs/cli/usage.md) for
the full CLI reference. See [docs/quickstart.md](docs/quickstart.md) for a
walkthrough with the bundled test data.

---

## Project Structure

Three top-level packages keep core logic separable from the UI runners:

```
MSpecViz/
├── mspecviz/                    # Core package — pure Python, no Streamlit imports
│   ├── search/                  # idXML parsing + peptide search
│   ├── spectra/                 # mzML streaming + Koina AI prediction
│   ├── annotation/              # Theoretical fragments + 7 comparison metrics
│   ├── dataframes/              # Polars LazyFrames for openms_insight
│   ├── plotting/                # Plotly + Matplotlib mirror/summary plots
│   ├── pipeline/                # Shared CLI/Streamlit orchestration
│   ├── export/                  # Shared export orchestrator + versioned dirs
│   └── common/                  # Ion colors, text helpers, parallel runner,
│                                # path disambiguation
├── streamlit_app/               # Streamlit UI (imports mspecviz.*)
│   ├── app.py                   # st.navigation entry point
│   ├── .streamlit/config.toml   # Streamlit server config (CWD-relative)
│   ├── config/
│   │   └── settings.json        # App-wide env config & enumerations
│   ├── utils/
│   │   ├── state/               # page_setup, workspaces, param I/O
│   │   └── page_state.py        # Session-state namespace helpers
│   └── views/
│       ├── setup.py             # Tab 1 — Setup
│       ├── search_load/         # Tab 2 — Search & Load
│       ├── explore/             # Tab 3 — Explore (3 sub-views)
│       ├── export/              # Tab 4 — Export
│       └── widgets/             # File browser, peptide input widget
├── cli/                         # Headless runner (imports mspecviz.*)
│   ├── main.py                  # Entry point
│   ├── pipeline.py              # Stage orchestration
│   ├── config.py                # CLIConfig dataclass + loaders
│   ├── config_builder.py        # 3-layer precedence merger + argparse
│   ├── exports.py               # DataFrame + plot export
│   └── config/                  # Example YAML configs
├── tests/                       # Unit + integration + GUI tests (pytest)
│   ├── cli/                     # CLI integration tests
│   ├── mspecviz/                # Core package tests (grouped by subpackage)
│   └── streamlit/               # GUI + export-pipeline tests
└── docs/                        # Reading guides (this README points into them)
```

**Rules the codebase follows:**

- `mspecviz/` must never import from `streamlit_app/` or `cli/` — the core
  stays testable and reusable.
- Both `streamlit_app/` and `cli/` depend on `mspecviz/`; anything shared
  between the two goes into `mspecviz/` (usually `mspecviz/pipeline/` for
  orchestration or `mspecviz/common/` for pure helpers).
- No `.py` file exceeds 250 lines; no folder contains more than 10 scripts.
  When a module grows past either cap, split it into a subpackage with a
  re-exporting `__init__.py`.

A per-subpackage tour lives in [docs/architecture/core.md](docs/architecture/core.md).

---

## Configuration at a glance

MSpecViz reads configuration from a small set of sources. Full reference
with every key documented is in
[docs/settings/overview.md](docs/settings/overview.md); here is the map:

| File | Scope | What it controls |
|---|---|---|
| `mspecviz/common/defaults.py` (`DEFAULTS`) | Single source of truth for tunable defaults | Every parameter's fallback value (ion types, tolerances, export DPI, Koina settings, …) |
| `streamlit_app/config/settings.json` | App-wide environmental config | Koina model catalogue, allowed ion types, data directory, performance tuning, ion color palette |
| `workspaces/<name>/params.json` | One per workspace | Sparse record of user-edited widget values; missing keys fall back to `DEFAULTS` |
| `cli/config/*.yaml` | One per CLI run | YAML config combining input files, peptides, annotation, export, and Koina sections |
| Command-line flags | Single invocation | Override any YAML field on a single run (`--workers`, `--mz-tolerance`, …) |

**Precedence for CLI runs:** `DEFAULTS` < YAML < CLI flags.

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

- `tests/mspecviz/` — unit + integration tests for the core (grouped into
  `loading/`, `annotation/`, `dataframes/`, `plotting/` subfolders)
- `tests/cli/` — CLI integration tests (invoke `python -m cli` as a
  subprocess)
- `tests/streamlit/` — GUI smoke tests (via `streamlit.testing.v1.AppTest`)
  plus export-pipeline unit tests

Marker-driven selection:

```bash
pytest tests/ -v -m "not integration"    # skip the slower end-to-end tests
pytest tests/mspecviz/annotation/ -v     # just the annotation subpackage
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
- Worker processes crash silently → re-run with `--workers 1` (CLI) or set
  the workers slider to 1 in the Setup tab to surface the real error.
