# Architecture Overview

How MSpecViz is layered, why the split exists, and what stays separable.
Detailed tours live in [core.md](core.md) and [frontends.md](frontends.md).

---

## Three-layer layout

```
MSpecViz/
├── core/              # Layer 1 — pure Python core, no UI imports
├── streamlit_app/     # Layer 2a — Streamlit UI, imports layer 1
├── cli/               # Layer 2b — headless runner, imports layer 1
└── tests/             # Layer 3 — pytest, imports all three
```

- **Layer 1** never imports from layer 2. Any helper that both the UI and
  the CLI need lives in layer 1 — orchestration in `core/pipeline/`,
  pure helpers in `core/helpers.py`.
- **Layers 2a and 2b are siblings** and never import from each other. They
  share code only through layer 1.
- **Layer 3** can import from everywhere. `tests/conftest.py` lives at the
  `tests/` root so every subfolder inherits the fixtures.

---

## Pipeline flow

A single run (CLI or Streamlit) flows through the modules in this order:

```
  config        cli/main.py:_load_config | streamlit_app/utils/state.py
    │           (YAML / settings.json / widgets)
    ▼
  lookup        core.peptide_lookup.load_idxmls → lookup_peptides
    │           → peptidoform_psms: dict[peptidoform → PSMList]
    ▼
  load+annotate core.pipeline.run_streaming_annotation
    │           (load_spectra_streaming + process_query_result:
    │            annotate, metrics, optional Koina fallback)
    │           → annotated_queries, reference_store, comparison_rows
    ▼
  dataframes    core.dataframes.build_psm_table + build_visualization_frames
    │           → psm / peaks / sequence / heatmap (Polars LazyFrames)
    ▼
  display/export Streamlit views/explore + views/export; both export
                paths delegate to core.export.run_export
                (versioned dir + parquets + plots + summaries).
```

The same `core.pipeline` helpers run in both entry points — only the parameter
source (YAML vs. session state) and output sink (disk vs. session state + disk)
differ.

---

## House rules

Checked by reviewers and by the repository's line/folder caps:

1. **≤250 lines per `.py`** — split into a subpackage when you cross it.
2. **≤10 scripts per folder** — subgroup by topic; shared helpers stay at the
   parent level.
3. **`core/` imports no UI code** — to cross that line, move the shared logic
   into `core/pipeline/` and have both callers pass data in.
4. **Re-export at the subpackage boundary** — callers import from
   `core.<subpackage>`, never a specific module inside.
5. **Module-level worker functions** — anything handed to a process pool must
   be picklable (no closures, no nested functions).
6. **Session-state prefix** — always `get_state_key(key)`; never write a bare
   string into `st.session_state`.

---

## See also

- [core.md](core.md) — module tour of `core/`.
- [frontends.md](frontends.md) — `streamlit_app/`, `cli/`, test layout.
- [../settings/overview.md](../settings/overview.md) — how configuration flows
  into this architecture.
