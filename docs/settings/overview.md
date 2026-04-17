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

## Key naming

Every setting has one name that is identical across all layers — the
`DEFAULTS` field, the `CLIConfig` attribute, the YAML key, the argparse
`dest`, the workspace `params.json` key, and the Streamlit session-state
key. No translation, no aliases.

Examples: `ion_types`, `mz_tolerance`, `add_neutral_losses`,
`export_directory`, `koina_enabled`, `koina_collision_energy`,
`figure_size`, `plot_unpaired_queries`, `summary_plot_types`.

See [cli.md](cli.md) for the full YAML schema and CLI flag table, and
`mspecviz/common/defaults.py` for the field list with default values.
Missing keys always fall back to `DEFAULTS.<field>`.
