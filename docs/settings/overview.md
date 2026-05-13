# Settings Overview

MSpecViz reads configuration from a small, explicit set of files. This
page shows *where* each source lives and *how* they are merged. For the
key-by-key reference, follow the links in each section.

---

## Where settings come from

| File / source | Scope | Who reads it | Detailed reference |
|---|---|---|---|
| `core/common/defaults.py` (`DEFAULTS`) | Tunable parameter defaults, single source of truth | Both front ends | (this file, see below) |
| `streamlit_app/config/settings.json` | App-wide environmental config & enumerations | Streamlit | [app.md](app.md) |
| `workspaces/<name>/params.json` | One workspace's live parameters | Streamlit only (written by `page_setup()`) | [workspace.md](workspace.md) |
| `cli/config/*.yaml` | One CLI invocation | CLI only | [cli.md](cli.md) |
| Command-line flags | One CLI invocation | CLI only | [cli.md](cli.md) |

All paths are relative to the **project root**. The Streamlit front end
`cd`s into `streamlit_app/` so it picks up `.streamlit/config.toml`, but
`workspaces/`, `exports/`, and config JSON still resolve against the
project root.

`DEFAULTS` is a frozen dataclass (`core.common.defaults.Defaults`).
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
DEFAULTS  <  YAML config
```

`cli/main.py:_load_config` builds a single dict from `DEFAULTS` (plus a
small `_CLI_EXTRAS` for CLI-only inputs), overlays the YAML, and wraps
the result in a `SimpleNamespace`. There are no command-line flag
overrides — per-run tweaks are made by editing the YAML.

---

## Naming convention

All keys are `snake_case` everywhere — Python attributes, JSON keys,
YAML keys, Streamlit session-state keys. Do not introduce kebab-case or
camelCase variants.

---

## Key naming

Every setting has one name that is identical across all layers — the
`DEFAULTS` field, the YAML key, the workspace `params.json` key, and the
Streamlit session-state key. No translation, no aliases.

Examples: `ion_types`, `mz_tolerance`, `add_neutral_losses`,
`export_directory`, `koina_enabled`, `koina_collision_energy`,
`figure_size`, `plot_unpaired_queries`, `summary_plot_types`.

See [cli.md](cli.md) for the full YAML schema, and
`core/common/defaults.py` for the field list with default values.
Missing keys always fall back to `DEFAULTS.<field>`.
