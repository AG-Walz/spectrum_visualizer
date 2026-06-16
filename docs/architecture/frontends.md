# Front-ends and Tests

The two UI runners that sit on top of `core/`, plus the test tree.

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
│   ├── progress_tracker.py     # Streamlit progress-bar adapter for run_tasks callbacks
│   └── state.py                # STATE_PREFIX, get_state_key, get_param, load_params,
│                               # save_params, workspace mgmt, page_setup, sidebar
└── views/
    ├── setup_page.py           # Tab 1 — file selection + parameters
    ├── load_annotate/
    │   ├── load_annotate_page.py    # Streamlit entry (render + fragments)
    │   └── load_annotate_actions.py # run_peptide_lookup, run_load_and_annotate
    ├── explore/
    │   ├── explore_page.py     # render + comparison_metrics_fragment
    │   ├── explorer.py         # openms_insight components wrapper
    │   └── mirror_plots.py     # mirror_plots_fragment
    ├── export/
    │   ├── export_page.py      # export_controls + history + zip download
    │   ├── export_actions.py   # _run_export (DataFrames + plots + settings)
    │   └── gallery.py          # _render_plot_gallery
    └── widgets/
        ├── file_browser.py     # Directory tree file picker
        └── peptide_input.py    # Peptide text-area + preview button
```

**Session-state convention.** Every key the app writes is prefixed with
`peptide_mirror.` (`STATE_PREFIX` in `utils/state.py` — a UI-layer concern,
not in core). Always resolve it through `get_state_key(key)`. Full key list:
[../settings/workspace.md#session-state-keys](../settings/workspace.md#session-state-keys).

**Page entry points.** Each `views/*/*_page.py` (and `setup_page.py`) exports a
`render()` called from `app.py`'s `st.navigation([st.Page(...)])` list.
Everything else in the folder is a support module (`*_actions.py` for
orchestration, `mirror_plots.py` / `gallery.py` for view-only pieces).

**Workspaces.** `utils/state.py` holds the workspace + sidebar helpers
(`save_params`, `is_safe_workspace_name`). Each workspace is a directory under
`workspaces/` with a sparse `params.json` of user-edited widget values; missing
keys fall back to `DEFAULTS`. Active workspace: `st.session_state["workspace"]`.

**Progress wiring.** `utils/progress_tracker.py` adapts the `(current, total,
message)` triples from `core.helpers` into `st.progress` calls — kept in the UI
layer so core stays framework-agnostic.

---

## CLI runner: `cli/`

```
cli/
├── __main__.py         # Enables `python -m cli <config.yaml>`
├── main.py             # YAML load + pipeline (lookup → load → export), tqdm progress
└── config/             # Example YAML configs (config_example.yaml + real-run configs)
```

The CLI is YAML-only: invoke `python -m cli <path/to/config.yaml>`. The
single positional argument is the config path; there are no flag overrides
— copy the YAML and edit per-run. Set `MSPECVIZ_LOG_LEVEL=ERROR` for quiet
runs.

`_load_config` merges `DEFAULTS` with the CLI-only `_CLI_EXTRAS`
(`idxml_files`, `mzml_files`, `peptide_source`, `peptides`,
`reference_substring`, `no_plots`) and then the YAML overrides into a
`SimpleNamespace`; `_expand_globs` and `_resolve_peptide_source` resolve the
file inputs. The merge order and the "add a parameter = one-field edit on
`DEFAULTS`" rule are documented in
[../settings/overview.md](../settings/overview.md).

Full reference: [../cli/usage.md](../cli/usage.md).

---

## Tests: `tests/`

```
tests/
├── conftest.py                       # fixtures (pytest walks up, subfolders inherit)
├── cli/
│   └── test_cli_integration.py       # subprocess-driven CLI integration tests
├── streamlit/
│   ├── test_gui.py                   # every view launches without exception (AppTest)
│   ├── test_export_pipeline.py
│   └── test_mirror_picker.py
└── core/                             # tests mirror the core layout
    ├── test_helpers.py
    ├── test_peak_annotation.py
    ├── test_comparison_metrics_unit.py
    ├── peptide_lookup/test_peptide_lookup.py
    ├── spectra/
    │   ├── test_koina_prediction.py
    │   └── test_spectrum_loader_integration.py
    ├── dataframes/
    │   ├── test_insight_frames_builders.py
    │   └── test_insight_frames_schemas.py
    ├── export/test_export_directories.py
    ├── pipeline/
    │   ├── test_load_annotate.py
    │   └── test_pipeline_characterization.py
    └── plotting/
        ├── test_mirror_plot.py
        ├── test_plot_worker.py
        └── test_summary_plots.py
```

Shared fixtures (defined in `tests/conftest.py`): `idxml_path`, `mzml_paths`,
`target_peptides`, `app_settings`, `spectrum_frame`, `make_peptide_match`.
Subfolder tests inherit them automatically; pytest walks the directory tree
upwards looking for `conftest.py`.

Marker `@pytest.mark.integration` is used for tests that touch real files
under `/mnt/volume/nacho/data/`.

Selective runs:

```bash
pytest tests/ -m "not integration"   # skip slow end-to-end tests
pytest tests/core/ -k "annotation or metrics"   # annotation + metrics units
pytest tests/cli/ -v                 # CLI integration
pytest tests/streamlit/ -v           # GUI smoke + export unit tests
```
