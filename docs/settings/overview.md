# Settings Overview

MSpecViz reads configuration from a small, explicit set of files. This
page shows *where* each source lives and *how* they are merged. For the
key-by-key reference, follow the links in each section.

---

## Where settings come from

| File / source | Scope | Who reads it | Detailed reference |
|---|---|---|---|
| `mspecviz/common/defaults.py` (`DEFAULTS`) | Tunable parameter defaults, single source of truth | Both front ends | (this file, see below) |
| `streamlit_app/config/settings.json` | App-wide environmental config & enumerations | Streamlit | [app.md](app.md) |
| `workspaces/<name>/params.json` | One workspace's live parameters | Streamlit only (written by `page_setup()`) | [workspace.md](workspace.md) |
| `cli/config/*.yaml` | One CLI invocation | CLI only | [cli.md](cli.md) |
| Command-line flags | One CLI invocation | CLI only | [cli.md](cli.md) |

All paths are relative to the **project root**. The Streamlit front end
`cd`s into `streamlit_app/` so it picks up `.streamlit/config.toml`, but
`workspaces/`, `exports/`, and config JSON still resolve against the
project root.

`DEFAULTS` is a frozen dataclass (`mspecviz.common.defaults.Defaults`).
Every tunable parameter that has a default value is defined there once;
the other layers only override. Change a default by editing that file.

---

## Precedence and merging

### Streamlit run

```
DEFAULTS  <  workspace params.json  <  live widget state
```

`settings.json` is loaded into `st.session_state.settings` at startup
and supplies enumerations (available Koina models, available ion types)
and environmental flags (workspace directory, cache sizing). It does
**not** contain tunable-parameter defaults — those come from `DEFAULTS`.

When a workspace is activated, its `params.json` is read directly.
Missing keys fall back to `DEFAULTS.<key>` at each call site. Any widget
change calls `save_params()` to write back immediately.

### CLI run

```
DEFAULTS  <  YAML config  <  argparse flags
```

`config_builder.build_config()` instantiates `CLIConfig` (whose field
defaults come from `DEFAULTS`) and applies the YAML and flag layers in
order. The mapping tables live in `cli/config_builder.py`
(`_YAML_OVERRIDES`, `_CLI_OVERRIDES`).

---

## Naming convention

All keys are `snake_case` everywhere — Python attributes, JSON keys,
YAML keys, Streamlit session-state keys, CLI flag destinations. Do not
introduce kebab-case or camelCase variants. CLI flags use kebab-case at
the shell (`--add-neutral-losses`) but argparse maps them to `snake_case`
destinations internally.

---

## Quick lookup: one setting, two places

The same conceptual setting appears in the workspace `params.json` and
the CLI YAML under slightly different shapes. The underlying default for
both comes from `DEFAULTS`.

| Concept | `DEFAULTS` field | Workspace `params.json` key | CLI YAML |
|---|---|---|---|
| Fragment ion types | `ion_types` | `ion_types` | `annotation.ion_types` |
| m/z tolerance | `mz_tolerance` | `mz_tolerance` | `annotation.mz_tolerance` |
| Neutral losses | `add_neutral_losses` | `add_neutral_losses` | `annotation.add_neutral_losses` |
| First prefix ion | `add_first_prefix_ion` | `add_first_prefix_ion` | `annotation.add_first_prefix_ion` |
| Remove precursor | `remove_precursor_peaks` | `remove_precursor_peaks` | `annotation.remove_precursor_peaks` |
| Intensity mode | `intensity_mode` | `intensity_mode` | `annotation.intensity_mode` |
| Workers | `workers` | `workers` | `performance.workers` |
| Export base dir | `export_directory` | `export_directory` | `export.directory` |
| Export formats | `export_formats` | *(not persisted)* | `export.formats` |
| Export DPI | `export_dpi` | *(not persisted)* | `export.export_dpi` |
| Figure size | `figure_size` | *(not persisted)* | `export.figure_size` |
| Plot unpaired queries | `plot_unpaired_queries` | *(not persisted)* | `export.plot_unpaired_queries` |
| Koina enabled | `koina_enabled` | `koina_enabled` | `koina.enabled` |
| Koina model | `koina_model` | `koina_model` | `koina.model` |
| Koina collision energy | `koina_collision_energy` | `koina_collision_energy` | `koina.collision_energy` |
| Koina instrument type | `koina_instrument_type` | `koina_instrument_type` | `koina.instrument_type` |
| Summary metric | `summary_metric` | *(not persisted)* | `summary.metric` |
| Summary plot types | `summary_plot_types` | *(not persisted)* | `summary.plot_types` |

The CLI YAML → `CLIConfig` mapping is encoded in
`cli/config_builder.py::_YAML_OVERRIDES`. The Streamlit workspace schema
is implicit: any key written by a widget is persisted; unknown keys are
preserved on reload. Missing keys always fall back to
`DEFAULTS.<field>`.
