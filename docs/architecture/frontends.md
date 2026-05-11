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
├── __main__.py         # Enables `python -m cli <config.yaml>`
├── main.py             # YAML load + pipeline (search → load → export), tqdm progress
└── config/             # Example YAML configs
```

The CLI is YAML-only: invoke `python -m cli <path/to/config.yaml>`. No
command-line flag overrides — copy the YAML and edit per-run.

**Two-layer config merge** (see `cli/main.py:_load_config`):

1. `mspecviz.common.defaults.DEFAULTS` + a small `_CLI_EXTRAS` dict
   (CLI-only fields: `idxml_files`, `mzml_files`, `peptide_source`,
   `reference_substring`, `no_plots`).
2. YAML file — per-run overrides (last-write-wins via dict-literal
   unpacking into a `SimpleNamespace`).

Adding a new parameter is a one-field edit on `DEFAULTS` (or
`_CLI_EXTRAS` for CLI-only inputs) — the YAML key inherits the name
automatically.

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
    ├── test_path_utils.py    # common/path_utils.py
    ├── test_export_dirs.py   # export/directories.py
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
`target_peptides`, `app_settings`, `spectrum_frame`, `make_peptide_match`.
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
