# Settings Overview

MSpecViz reads configuration from up to four places depending on the run
mode. This page shows *where* each source lives and *how* they are merged.
For the key-by-key reference, follow the links in each section.

---

## Where settings come from

| File / source | Scope | Who reads it | Detailed reference |
|---|---|---|---|
| `streamlit_app/config/settings.json` | App-wide, rarely changed | Both front ends (CLI pulls layer-1 defaults here too) | [app.md](app.md) |
| `streamlit_app/config/default-parameters.json` | Initial values for a new workspace | Streamlit only | [workspace.md](workspace.md) |
| `workspaces/<name>/params.json` | One workspace's live parameters | Streamlit only (written by `page_setup()`) | [workspace.md](workspace.md) |
| `cli/config/*.yaml` | One CLI invocation | CLI only | [cli.md](cli.md) |
| Command-line flags | One CLI invocation | CLI only | [cli.md](cli.md) |

All paths are relative to the **project root**. The Streamlit front end
`cd`s into `streamlit_app/` so it picks up `.streamlit/config.toml`, but
`workspaces/`, `exports/`, and config JSON still resolve against the
project root.

---

## Precedence and merging

### Streamlit run

```
settings.json  <  default-parameters.json  <  workspace params.json  <  live widget state
```

`settings.json` values are loaded into `st.session_state.settings` at
startup. When a workspace is activated, `default-parameters.json` is
copied into `workspaces/<name>/params.json` (if it doesn't already exist);
`params.json` then becomes the source of truth for user-editable
parameters. Any widget change calls `save_params()` to write back
immediately.

### CLI run

```
settings.json  <  YAML config  <  argparse flags
```

`config_builder.build_config()` merges the three layers in order. Each
layer is a small table in `cli/config_builder.py`. Missing values fall
back to the dataclass defaults on `CLIConfig`.

---

## Quick lookup: one setting, three places

The same conceptual setting surfaces differently in each configuration
system. Note the different key-casing conventions:

- `settings.json` — nested dotted sections, `snake_case` leaves.
- `default-parameters.json` — flat `hyphen-case` (also the Streamlit
  widget keys).
- CLI YAML — nested dotted sections, `snake_case` leaves (same shape as
  `settings.json`).

| Concept | `settings.json` | `default-parameters.json` | CLI YAML |
|---|---|---|---|
| Fragment ion types | `ion_types.default` | `ion-types` | `annotation.ion_types` |
| m/z tolerance | `mz_tolerance.default` | `mz-tolerance` | `annotation.mz_tolerance` |
| Neutral losses | `annotation.add_losses` | `add-losses` | `annotation.add_losses` |
| First prefix ion | `annotation.add_first_prefix_ion` | `add-first-prefix-ion` | `annotation.add_first_prefix_ion` |
| Remove precursor | `annotation.remove_precursor_peaks` | `remove-precursor-peaks` | `annotation.remove_precursor_peaks` |
| Workers | `performance.workers` | `workers` | `performance.workers` |
| Export base dir | `export.default_directory` | `export-dir` | `export.directory` |
| Export formats | `export.formats` | *(not in file; live widget state)* | `export.formats` |
| DPI | `export.dpi` | *(not in file; live widget state)* | `export.dpi` |
| Figsize | `export.figsize` | *(not in file; live widget state)* | `export.figsize` |
| Koina enabled | `koina.enabled` | `koina-enabled` | `koina.enabled` |
| Koina model | `koina.model` | `koina-model` | `koina.model` |
| Koina CE | `koina.collision_energy` | `koina-ce` | `koina.collision_energy` |
| Koina instrument | `koina.instrument_type` | `koina-instrument` | `koina.instrument_type` |

The three-way mapping is encoded in `cli/config_builder.py`
(`_SETTINGS_DEFAULTS`, `_YAML_OVERRIDES`) and
`streamlit_app/utils/state/params.py`.
