# Workspace & Session-State Settings

Two sources shape the Streamlit per-workspace state:

- `workspaces/<name>/params.json` — one per workspace, persisted widget state
- `peptide_mirror.*` keys in `st.session_state` — live widget state

Defaults come from `core.defaults.DEFAULTS` (see
[overview.md](overview.md)). There is no template file: a new workspace
starts with an empty `params.json` and each widget resolves its initial
value via `get_param(params, key, DEFAULTS.<key>)`.

---

## `workspaces/<name>/params.json`

Written by `utils/state.py::save_params()` whenever a
widget value changes. Keys are `snake_case` and match the widget keys
used in session state (stored as `peptide_mirror.<key>` — the prefix is
stripped for on-disk persistence). The schema is sparse: only keys the
user has touched get persisted, and any key missing from the file falls
back to `DEFAULTS.<key>` at each call site.

Most persisted keys are tunable parameters whose names, types, and
default values live in `core/defaults.py` (`ion_types`, `mz_tolerance`,
`add_neutral_losses`, `intensity_mode`, the `koina_*` group, the
`summary_*` group, `export_directory`, `project_name`, `workers`,
`image_format`, …). That dataclass is the single source of truth — this
doc does not re-list them. The only `params.json` keys **not** backed by
`DEFAULTS` are the Setup-tab inputs, which have no other home:

| Key | Type | Fallback | Description |
|---|---|---|---|
| `selected_idxmls` | list[string] | `[]` | Absolute paths of idXML files queued in the Setup tab. |
| `selected_mzmls` | list[string] | `[]` | Absolute paths of mzML files queued in the Setup tab. |
| `idxml_browser_dir` | string | `""` | Last directory browsed for idXML selection. |
| `mzml_browser_dir` | string | `""` | Same, for the mzML browser. |
| `target_peptides` | string | `""` | Raw contents of the peptide-input textarea (newline-separated). Parsed by `parse_peptide_list` at lookup time. |
| `reference_substring` | string | `""` | Filename substring identifying reference mzMLs. Empty ⇒ single-spectrum mode, no mirror plots. |

**Editing defaults.** Change `core/defaults.py`. This affects
every workspace that has not explicitly overridden the key (and every
CLI run that does not override it). Existing workspaces that did
override a key keep their value in `params.json`.

File layout under the project root:

```
workspaces/
├── default/
│   └── params.json
├── alice_lumos_hcd/
│   └── params.json
└── alice_timstof/
    └── params.json
```

Workspace names are validated by `is_safe_workspace_name()` in
`streamlit_app/utils/state.py` — no `..`, no absolute
paths, no shell metacharacters.

---

## Session-state keys

Every Streamlit widget writes into `st.session_state` under the
`peptide_mirror.` namespace (`STATE_PREFIX` in
`streamlit_app/utils/state.py` — the prefix is a UI-layer
concern, not part of the core package). Use `get_state_key(key)` to
resolve the prefixed version.

The parameter keys above (the `params.json` set) also live in session
state while the app runs — Setup writes them, Load & Annotate and Export
read them. The keys that matter as a cross-view contract are the
**derived artifacts** that one stage produces and a later view consumes:

| Key | Written by | Read by | Description |
|---|---|---|---|
| `peptidoform_psms` | Load & Annotate (stage 1) | Load & Annotate (stage 2), Export | Lookup results — the `Dict[peptidoform, PSMList]` returned by `lookup_peptides`. |
| `psm_frame`, `psm_lookup`, `peaks_frame`, `sequence_frame`, `heatmap_frame` | Load & Annotate (stage 2) | Explore, Export | Polars LazyFrames for openms_insight. |
| `annotated_queries`, `reference_store`, `comparison_frame` | Load & Annotate (stage 2) | Explore, Export | Annotated spectra + metrics. |
| `data_ready` | Load & Annotate (stage 2) | Explore, Export | `True` when all four LazyFrames are populated — gates the downstream views. |
| `cache_version` | Load & Annotate (stage 2) | Cache validator | Bumped whenever a fresh load completes; used to detect re-runs and invalidate cached frames. |

---

## See also

- [app.md](app.md) — the app-wide `settings.json` file.
- [cli.md](cli.md) — CLI-only YAML reference (the CLI is YAML-only — no
  flag overrides).
- [../architecture/frontends.md](../architecture/frontends.md) — where
  these session-state keys get read and written in the views.
