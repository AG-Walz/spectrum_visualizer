# App-wide settings — `settings.json`

Path: `streamlit_app/config/settings.json`

App-wide environmental configuration and enumerations. Read once at
startup into `st.session_state.settings` (`streamlit_app/app.py`).
Top-level keys are read directly (`settings["app_name"]`); nested keys
go through `_get_nested_value(settings, "section", "subkey",
default=...)` in `streamlit_app/views/setup_page.py`, so a missing key
falls back to a hard-coded default.

`settings.json` does **not** own tunable-parameter defaults. Those live
in `core/defaults.py` (`DEFAULTS`) — see [overview.md](overview.md).
This file is Streamlit-only; the CLI reads `DEFAULTS` directly and never
opens `settings.json`.

---

## Keys

The file has exactly five top-level keys:

| Key | Type | Description |
|---|---|---|
| `app_name` | string | Shown in the Streamlit page title (`st.set_page_config`). |
| `enable_workspaces` | bool | If `false`, the sidebar workspace switcher is hidden and all state lives in a single default workspace. |
| `workspaces_dir` | string | Directory (relative to project root) where workspace folders are created. |
| `data_directories.default` | string | Starting directory for the Setup tab's file browser. Point at your lab's shared storage root; symlinks are followed. (`._comment` is a non-functional note.) |
| `ion_types.available` | list[string] | Full set of ion series the UI offers. Values outside this list are ignored by the annotation step. The pre-selected subset comes from `DEFAULTS.ion_types`. |

---

## Facts that live elsewhere

- **Koina model catalogue** → `core/defaults.py:KOINA_MODELS` (imported
  directly by `streamlit_app/views/setup_page.py`).
- **Ion colour palette** → `core/defaults.py:ION_COLORS` (used by
  `core/plotting/mirror.py`).
- **Tunable-parameter defaults** (mz tolerance, workers, cache, etc.)
  → `core/defaults.py:DEFAULTS`.
