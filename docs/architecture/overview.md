# Architecture Overview

How MSpecViz is layered, why the split exists, and what stays separable.
Detailed tours live in [core.md](core.md) and [frontends.md](frontends.md).

---

## Three-layer layout

```
MSpecViz/
├── mspecviz/          # Layer 1 — pure Python core, no UI imports
├── streamlit_app/     # Layer 2a — Streamlit UI, imports layer 1
├── cli/               # Layer 2b — headless runner, imports layer 1
└── tests/             # Layer 3 — pytest, imports all three
```

- **Layer 1** never imports from layer 2. Any helper that both the UI and
  the CLI need lives in layer 1 — orchestration in `mspecviz/pipeline/`,
  pure helpers in `mspecviz/common/`.
- **Layers 2a and 2b are siblings** and never import from each other. They
  share code only through layer 1.
- **Layer 3** can import from everywhere. `tests/conftest.py` lives at the
  `tests/` root so every subfolder inherits the fixtures.

---

## Pipeline flow

A single run (CLI or Streamlit) flows through the subpackages in this order:

```
  YAML / settings.json / widgets
          │
          ▼
  [ config layer ]            cli/main.py:_load_config, utils/state/params.py
          │
          ▼
  [ search ]                  mspecviz.search.search_peptides_in_idxmls
          │   → variant_groups: dict[variant_key → VariantGroup]
          │
          ▼
  [ load + annotate ]         mspecviz.spectra.load_spectra_streaming
          │                   mspecviz.pipeline.process_query_result
          │                   (invokes annotation + metrics + Koina fallback)
          │   → annotated_queries, reference_store, comparison_frame
          │
          ▼
  [ dataframes ]              mspecviz.dataframes.build_psm_table
          │                   build_peaks_dataframe / sequence / heatmap
          │   → psm_frame, peaks_frame, sequence_frame, heatmap_frame (Polars LazyFrames)
          │
          ▼
  [ display / export ]        Streamlit: views/explore + views/export
                              Both export paths delegate to
                              mspecviz.export.run_export
                              (versioned dir + parquets + plots + summaries).
```

The same `mspecviz.pipeline` helpers run in both entry points — the only
difference is where parameters come from (YAML vs. session state) and
where outputs go (stdout/disk vs. session state + disk).

---

## House rules

Checked by reviewers and by the repository's line/folder caps:

1. **≤250 lines per `.py`** — split into a subpackage when you cross it.
2. **≤10 scripts per folder** — subgroup by topic (`search/models.py`,
   `search/workers.py`, …). Shared helpers stay at the parent level.
3. **`mspecviz/` imports no UI code** — if you need to cross that line,
   move the shared logic into `mspecviz/pipeline/` and have both callers
   pass data in.
4. **Re-export at the subpackage boundary** — external callers import from
   `mspecviz.<subpackage>`, never from a specific module inside.
5. **Module-level worker functions** — any function handed to
   `ProcessPoolExecutor` must be picklable, which means no closures and no
   nested functions. Put them at module top level.
6. **Session-state prefix** — always `get_state_key(key)`. Never write a
   bare string into `st.session_state`.

---

## See also

- [core.md](core.md) — per-subpackage tour of `mspecviz/`.
- [frontends.md](frontends.md) — `streamlit_app/`, `cli/`, and test layout.
- [../settings/overview.md](../settings/overview.md) — how configuration
  flows into this layered architecture.
