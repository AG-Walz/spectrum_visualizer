# Workspace & Session-State Settings

Two closely related files and the live Streamlit session state:

- `streamlit_app/config/default-parameters.json` — template
- `workspaces/<name>/params.json` — one per workspace, materialized from
  the template above
- `peptide_mirror.*` keys in `st.session_state` — live widget state

---

## `streamlit_app/config/default-parameters.json`

Seeds every new workspace's `params.json`. These are the values a user
typically changes at the start of a new analysis — input files, target
peptides, tolerances, Koina flags. Keys are the same ones used as widget
keys throughout the app and stored in session state as
`peptide_mirror.<key>`.

| Key | Type | Default | Description |
|---|---|---|---|
| `selected-idxmls` | list[string] | `[]` | Absolute paths of idXML files queued in the Setup tab. |
| `selected-mzmls` | list[string] | `[]` | Absolute paths of mzML files queued in the Setup tab. |
| `idxml-browser-dir` | string | `""` | Last directory browsed for idXML selection. Persisted so the Setup tab opens where you left off. |
| `mzml-browser-dir` | string | `""` | Same, for the mzML browser. |
| `target-peptides` | string | `""` | Raw contents of the peptide-input textarea (newline-separated). Parsed by `parse_peptide_list` at search time. |
| `reference-substring` | string | `""` | Filename substring identifying reference mzMLs (e.g. `Synpep`). Empty ⇒ single-spectrum mode, no mirror plots. |
| `ion-types` | list[string] | `["b", "y"]` | Active fragment ion series. Must be a subset of `settings.json > ion_types.available`. |
| `mz-tolerance` | float | `0.02` | m/z matching tolerance in Da. |
| `add-losses` | bool | `false` | Include neutral-loss ions in annotation. |
| `add-first-prefix-ion` | bool | `true` | Include a1/b1/c1. |
| `remove-precursor-peaks` | bool | `false` | Strip precursor peaks before annotation. |
| `export-dir` | string | `""` | Base export directory. Empty ⇒ `settings.json > export.default_directory` is used. |
| `workers` | int | `4` | Parallelism override for this workspace. |
| `intensity-mode` | string | `"Relative"` | `"Relative"` (normalize to base peak = 1.0) or `"Raw"` (preserve original intensities). Affects both plots and metrics. |
| `koina-enabled` | bool | `false` | Per-workspace Koina flag. Note: `settings.json > koina.enabled` only controls whether the UI *exposes* Koina at all; this key controls whether it is *used* for the current workspace. |
| `koina-model` | string | `"ms2pip_timsTOF2024"` | Koina model for this workspace. |
| `koina-ce` | int | `35` | Collision energy. |
| `koina-instrument` | string | `"TIMSTOF"` | Instrument type for Koina. |
| `image-format` | string | `"png"` | Default Streamlit download format for inline plot previews. |

**Editing defaults.** Modifying `default-parameters.json` only affects
**new** workspaces. Existing workspaces keep their `params.json`. To
propagate a change, either delete the workspace's `params.json` or edit
it directly.

---

## `workspaces/<name>/params.json`

Written by `utils/state/params.py`'s `save_params()` whenever a widget
value changes. Same schema as `default-parameters.json`. You can edit it
by hand for reproducible runs; the Streamlit app reloads it on next
launch. Versioned export directories reference a fingerprint of this file
for reuse detection.

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
`streamlit_app/utils/state/workspace.py` — no `..`, no absolute paths, no
shell metacharacters.

---

## Session-state keys

Every Streamlit widget writes into `st.session_state` under the
`peptide_mirror.` namespace (`STATE_PREFIX` in
`mspecviz/common/constants.py`). The keys below are the public contract
between views — written by one view, read by another.

| Key | Written by | Read by | Description |
|---|---|---|---|
| `selected_idxmls`, `selected_mzmls`, `target_peptides`, `reference_substring` | Setup | Search & Load | Input definitions. |
| `ion_types`, `mz_tolerance`, `add_losses`, `add_first_prefix_ion`, `remove_precursor_peaks`, `intensity_mode` | Setup | Search & Load, Export | Annotation parameters. |
| `workers`, `koina_enabled`, `koina_model`, `koina_ce`, `koina_instrument` | Setup | Search & Load, Export | Runtime / AI-prediction parameters. |
| `variant_groups` | Search & Load (stage 1) | Search & Load (stage 2), Export | Search results — the `dict[variant_key → VariantGroup]` returned by `search_peptides_in_idxmls`. |
| `psm_lf`, `psm_lookup`, `peaks_lf`, `seq_lf`, `heatmap_lf` | Search & Load (stage 2) | Explore, Export | Polars LazyFrames for openms_insight. |
| `annotated_queries`, `reference_store`, `comparison_df` | Search & Load (stage 2) | Explore, Export | Annotated spectra + metrics. |
| `data_ready` | Search & Load (stage 2) | Explore, Export | `True` when all four LazyFrames are populated — gates the downstream views. |
| `cache_version` | `page_setup()` | Cache validator | Bumped whenever the active workspace changes or the user clicks "Invalidate cache". |

All keys are prefixed with `peptide_mirror.` when actually stored. Use
`get_state_key(key)` to resolve the prefixed version.

---

## See also

- [app.md](app.md) — the app-wide `settings.json` file.
- [cli.md](cli.md) — CLI-only YAML + flag reference.
- [../architecture/frontends.md](../architecture/frontends.md) — where
  these session-state keys get read and written in the views.
