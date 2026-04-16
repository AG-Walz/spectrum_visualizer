# Front-ends and Tests

The two UI runners that sit on top of `mspecviz/`, plus the test tree.

---

## Streamlit front end: `streamlit_app/`

```
streamlit_app/
├── app.py                      # st.navigation + st.Page entry point
├── .streamlit/config.toml      # server config (picked up from CWD)
├── config/
│   └── settings.json           # app-wide env config & enumerations
├── utils/
│   ├── page_state.py           # session-state namespace helpers
│   └── state/                  # page_setup, workspace mgmt, sidebar
│       ├── core.py             # get_state_key, _persisted_key, get_param
│       ├── workspace.py        # is_safe_workspace_name, switcher UI
│       └── params.py           # load_params, save_params, render_sidebar
└── views/
    ├── setup.py                # Tab 1 — file selection + parameters
    ├── search_load/
    │   ├── page.py             # Streamlit entry (render + fragments)
    │   └── pipeline.py         # run_peptide_search, run_search_and_load
    ├── explore/
    │   ├── page.py             # render + comparison_metrics_fragment
    │   ├── explorer.py         # openms_insight components wrapper
    │   └── plots.py            # mirror_plots_fragment
    ├── export/
    │   ├── page.py             # export_controls + history + zip download
    │   ├── pipeline.py         # _run_export (DataFrames + plots + settings)
    │   └── gallery.py          # _render_plot_gallery
    └── widgets/
        ├── file_browser.py     # Directory tree file picker
        └── peptide_input.py    # Peptide text-area + preview button
```

**Session-state convention.** Every key the app writes is prefixed with
`peptide_mirror.` (see `STATE_PREFIX` in `mspecviz/common/constants.py`).
Always go through `get_state_key(key)` — never hard-code the prefix. Full
list of public keys: [../settings/workspace.md#session-state-keys](../settings/workspace.md#session-state-keys).

**Page entry points.** Each `views/*/page.py` (or `setup.py`) exports a
`render()` function called from `app.py`'s `st.navigation([st.Page(...)])`
list. Everything else in the folder is a support module.

**Workspaces.** Managed by `utils/state/workspace.py`. Each workspace is a
directory under `workspaces/` with its own `params.json` — a sparse
record of user-edited widget values. Missing keys fall back to
`mspecviz.common.defaults.DEFAULTS`. The active workspace lives in
`st.session_state["workspace"]`.

---

## CLI runner: `cli/`

```
cli/
├── __main__.py         # Enables `python -m cli`
├── main.py             # Argparse + run() entry point
├── config.py           # CLIConfig dataclass + YAML/peptide loaders
├── config_builder.py   # 3-layer precedence merge + argparse surface
├── pipeline.py         # _validate_inputs, _run_search, _run_load, run_pipeline
├── exports.py          # export_dataframes, prepare_export_dir, export_plots
├── output.py           # TerminalProgress, log formatting
├── plots.py            # CLI-facing plot wrapper (calls mspecviz.plotting)
└── config/             # Example YAML configs
```

**Three-layer config merge** (see `cli/config_builder.py`):

1. `streamlit_app/config/settings.json` — app defaults
2. YAML config file — per-run overrides
3. argparse flags — single-invocation overrides (highest priority)

Each layer is a table in `config_builder.py` (`_SETTINGS_DEFAULTS`,
`_YAML_OVERRIDES`, `_CLI_OVERRIDES`). Adding a new parameter is a
one-row edit to each relevant table, not three scattered branches.

Full reference: [../cli/usage.md](../cli/usage.md).

---

## Tests: `tests/`

```
tests/
├── conftest.py               # fixtures (pytest walks up, subfolders inherit)
├── cli/                      # subprocess-driven CLI integration tests
├── streamlit/                # AppTest GUI smoke + export unit tests
│   ├── test_gui.py           # every view launches without exception
│   ├── test_export_pipeline.py
│   └── test_export_bugfixes.py
└── mspecviz/
    ├── test_path_utils.py    # parent-level shared helpers
    ├── test_export_dirs.py
    ├── loading/              # search + spectrum loader + Koina
    ├── annotation/           # peak annotation + spectrum metrics
    │   ├── test_peak_annotation_core.py
    │   ├── test_peak_annotation_integration.py
    │   ├── test_spectrum_metrics_unit.py
    │   └── test_spectrum_metrics_integration.py
    ├── dataframes/           # schema + builder tests
    └── plotting/             # mirror + summary + worker tests
```

Shared fixtures (defined in `tests/conftest.py`): `idxml_path`, `mzml_paths`,
`target_peptides`, `app_settings`, `spectrum_df`, `make_peptide_match`.
Subfolder tests inherit them automatically; pytest walks the directory tree
upwards looking for `conftest.py`.

Marker `@pytest.mark.integration` is used for tests that touch real files
under `/mnt/volume/nacho/data/`.

Selective runs:

```bash
pytest tests/ -m "not integration"  # skip slow end-to-end tests
pytest tests/mspecviz/annotation/    # just one subpackage
pytest tests/cli/ -v                 # CLI integration
```
