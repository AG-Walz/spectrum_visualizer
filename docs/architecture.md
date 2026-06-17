# Architecture

How MSpecViz is layered, why the split exists, and the rules that keep it
separable. This describes structure and rationale — for the current file list,
browse the code; for tunables, `core/defaults.py` is the source of truth.

---

## Three layers

- **`core/`** — pure Python, no UI imports. The pipeline, helpers, and every
  tunable default live here.
- **`streamlit_app/`** and **`cli/`** — sibling front ends. Both import `core/`;
  neither imports the other. Any logic they share moves into `core/`
  (orchestration into `core/pipeline/`, pure helpers into `core/helpers.py`).
- **`tests/`** — imports all three. `tests/conftest.py` sits at the root so every
  subfolder inherits its fixtures.

`core/` never importing UI code is the load-bearing rule: it keeps the pipeline
testable and reusable from both runners.

---

## Pipeline flow

One run (CLI or Streamlit) flows through the same `core` helpers in order:

1. **config** — `cli/main.py:_load_config` (YAML) or `streamlit_app/utils/state.py`
   (session state). Only the source differs.
2. **lookup** — `core.peptide_lookup`: `load_idxmls` → `lookup_peptides`, yielding
   `dict[peptidoform → PSMList]`.
3. **load + annotate** — `core.pipeline.run_streaming_annotation`: stream spectra,
   annotate, compute metrics, optional Koina fallback.
4. **dataframes** — `core.dataframes`: `build_psm_table` + `build_visualization_frames`
   → the four Polars LazyFrames openms_insight consumes (see [export.md → data files](export.md#data-files)).
5. **display / export** — Streamlit `views/` or `core.export.run_export` (both
   export paths call the same orchestrator, so outputs are byte-identical; see
   [export.md → directory layout and versioning](export.md#directory-layout-and-versioning)).

Only the parameter source (YAML vs. session state) and output sink (disk vs.
session state + disk) differ between the two front ends.

---

## `core/` at a glance

Flat top-level modules plus three subpackages. Each module's own docstring and
constants are the source of truth; below is just the map and the non-obvious
rationale.

**Top-level modules:** `defaults.py` (the `DEFAULTS` dataclass + `KOINA_MODELS`,
`INSTRUMENT_TYPES`, `ION_COLORS`), `helpers.py` (`run_tasks` — the single parallel
dispatcher wrapping `joblib`, plus text/path helpers), `peptide_lookup.py` (the
only place idXML is parsed — see below), `peak_annotation.py` (theoretical
fragments via `spectrum_utils.proforma`), `comparison_metrics.py`
(`compute_all_metrics` — see [metrics.md → the seven metrics](metrics.md#the-seven-metrics)),
`dataframes.py` (the four LazyFrame builders + schemas — see [export.md → data files](export.md#data-files)),
`export.py` (the shared, versioned export orchestrator — see
[export.md → directory layout and versioning](export.md#directory-layout-and-versioning)).

**`core/spectra/`** — mzML streaming + Koina prediction, both producing one
`LoadedSpectrum` so the annotation layer ignores origin. Native-ID normalisation
lives in `mzml_native_ids.py`: **timsTOF/Bruker `index=` IDs are 1-based** and the
loader auto-detects the offset.

**`core/plotting/`** — one mirror renderer (`mirror.py`, branching Matplotlib vs.
Plotly), the picklable parallel plot worker (`plot_worker.py`), and summary plots
(`summary.py`). Two rationale notes: pyopenms-viz' `plot_spectrum` was abandoned
because it scrambles ion colouring when ion types interleave; and only `dual_axis`
and `ridgeline` summary plots exist — the Explore-tab Heatmap is a separate
interactive widget, never exported (see [export.md → plots](export.md#plots)).

**`core/pipeline/`** — the glue the CLI and Streamlit share. Any feature that must
run identically from both front ends belongs here. `kwargs_builders.py` takes a
generic `getter(key, default)` so each front end supplies parameters without core
knowing its storage format.

`load_idxmls` parses each idXML once into an in-memory `PSMList`; `lookup_peptides`
and `count_psms_per_peptide` then scan that list — the same file is never read
twice.

---

## Front ends

### Streamlit — `streamlit_app/`

Imports only from `core/`. `app.py` builds the `st.navigation([st.Page(...)])` list
and nothing else. Each tab is a folder under `views/` (`setup_page.py`, plus
`load_annotate/`, `explore/`, `export/`, and shared `widgets/`); every page module
exports `render()`, and `*_actions.py` siblings hold orchestration.

`utils/state.py` is the hub: `STATE_PREFIX`/`get_state_key`, the param accessors
(`get_param`, `load_params`, `save_params`), workspace management, `page_setup`, and
the sidebar. `progress_tracker.py` adapts the `(current, total, message)` triples
from `core.helpers` into `st.progress` so core stays framework-agnostic.

Session-state and workspace details: [settings.md → workspaces and session state](settings.md#workspaces-and-session-state).

### CLI — `cli/`

`python -m cli <config.yaml>` (enabled by `__main__.py`) runs the 3-step pipeline
in `main.py` with tqdm progress. YAML-only: one positional config path, no flag
overrides. Full reference: [cli.md → invocation](cli.md#invocation).

### Tests — `tests/`

The tree mirrors `core/`. `tests/conftest.py` holds the shared fixtures
(`idxml_path`, `mzml_paths`, `target_peptides`, `app_settings`, `spectrum_frame`,
`make_peptide_match`) and pytest walks upward, so subfolders inherit them.
`tests/cli/` runs subprocess-driven CLI tests; `tests/streamlit/` runs an AppTest
GUI smoke test plus export/picker units. `@pytest.mark.integration` marks tests
touching real files under `/mnt/volume/nacho/data/`.

```bash
pytest tests/ -m "not integration"              # skip slow end-to-end tests
pytest tests/core/ -k "annotation or metrics"   # annotation + metrics units
```
