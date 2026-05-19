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

Written by `utils/state/state_params.py::save_params()` whenever a
widget value changes. Keys are `snake_case` and match the widget keys
used in session state (stored as `peptide_mirror.<key>` — the prefix is
stripped for on-disk persistence). The schema is sparse: only keys the
user has touched get persisted. Missing keys fall back to
`DEFAULTS.<key>` at each call site.

| Key | Type | Fallback | Description |
|---|---|---|---|
| `selected_idxmls` | list[string] | `[]` | Absolute paths of idXML files queued in the Setup tab. |
| `selected_mzmls` | list[string] | `[]` | Absolute paths of mzML files queued in the Setup tab. |
| `idxml_browser_dir` | string | `""` | Last directory browsed for idXML selection. |
| `mzml_browser_dir` | string | `""` | Same, for the mzML browser. |
| `target_peptides` | string | `""` | Raw contents of the peptide-input textarea (newline-separated). Parsed by `parse_peptide_list` at search time. |
| `reference_substring` | string | `""` | Filename substring identifying reference mzMLs. Empty ⇒ single-spectrum mode, no mirror plots. |
| `ion_types` | list[string] | `DEFAULTS.ion_types` | Active fragment ion series. Must be a subset of `settings.json > ion_types.available`. |
| `mz_tolerance` | float | `DEFAULTS.mz_tolerance` | m/z matching tolerance in Da. |
| `add_neutral_losses` | bool | `DEFAULTS.add_neutral_losses` | Include neutral-loss ions in annotation. |
| `export_directory` | string | `DEFAULTS.export_directory` (`"./exports"`) | Base export directory. |
| `project_name` | string | `DEFAULTS.project_name` | Prefix for the versioned export directory. |
| `workers` | int | `DEFAULTS.workers` | Parallelism override for this workspace. |
| `intensity_mode` | string | `DEFAULTS.intensity_mode` | `"Relative"` or `"Raw"`. Affects both plots and metrics. |
| `koina_enabled` | bool | `DEFAULTS.koina_enabled` | Per-workspace Koina flag. |
| `koina_model` | string | `DEFAULTS.koina_model` | Koina model for this workspace. |
| `koina_collision_energy` | int | `DEFAULTS.koina_collision_energy` | Collision energy. |
| `koina_instrument_type` | string | `DEFAULTS.koina_instrument_type` | Instrument type for Koina. |
| `image_format` | string | `DEFAULTS.image_format` | Streamlit single-image picker default (PNG/SVG/PDF/HTML). CLI uses `export_formats` (plural) instead. |
| `summary_metric` | string | `DEFAULTS.summary_metric` | Metric the summary plots are coloured/sorted by. |
| `summary_plot_types` | list[string] | `DEFAULTS.summary_plot_types` | Which summary plots to generate (`dual_axis`, `ridgeline`). |

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
`streamlit_app/utils/state/state_workspace.py` — no `..`, no absolute
paths, no shell metacharacters.

---

## Session-state keys

Every Streamlit widget writes into `st.session_state` under the
`peptide_mirror.` namespace (`STATE_PREFIX` in
`streamlit_app/utils/state/state_core.py` — the prefix is a UI-layer
concern, not part of the core package). The keys below are the public
contract between views — written by one view, read by another.

| Key | Written by | Read by | Description |
|---|---|---|---|
| `selected_idxmls`, `selected_mzmls`, `target_peptides`, `reference_substring` | Setup | Search & Load | Input definitions. |
| `ion_types`, `mz_tolerance`, `add_neutral_losses`, `intensity_mode` | Setup | Search & Load, Export | Annotation parameters. |
| `workers`, `koina_enabled`, `koina_model`, `koina_collision_energy`, `koina_instrument_type` | Setup | Search & Load, Export | Runtime / AI-prediction parameters. |
| `variant_groups` | Search & Load (stage 1) | Search & Load (stage 2), Export | Search results — the `Dict[variant_key, PSMList]` returned by `search_peptides_in_idxmls`. |
| `psm_frame`, `psm_lookup`, `peaks_frame`, `sequence_frame`, `heatmap_frame` | Search & Load (stage 2) | Explore, Export | Polars LazyFrames for openms_insight. |
| `annotated_queries`, `reference_store`, `comparison_frame` | Search & Load (stage 2) | Explore, Export | Annotated spectra + metrics. |
| `data_ready` | Search & Load (stage 2) | Explore, Export | `True` when all four LazyFrames are populated — gates the downstream views. |
| `cache_version` | Search & Load (stage 2) | Cache validator | Bumped whenever a fresh load completes; used to detect re-runs and invalidate cached frames. |

All keys are prefixed with `peptide_mirror.` when actually stored. Use
`get_state_key(key)` to resolve the prefixed version.

---

## See also

- [app.md](app.md) — the app-wide `settings.json` file.
- [cli.md](cli.md) — CLI-only YAML reference (the CLI is YAML-only — no
  flag overrides).
- [../architecture/frontends.md](../architecture/frontends.md) — where
  these session-state keys get read and written in the views.
