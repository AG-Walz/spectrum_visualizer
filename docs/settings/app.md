# App-wide settings — `settings.json`

Path: `streamlit_app/config/settings.json`

App-wide settings — changed rarely and centrally. All keys are read via
`get_nested(settings, "section", "subkey", default=...)` so missing keys
fall back to a hard-coded default.

Both the Streamlit app and the CLI read from this file. The CLI uses it
as the first layer in its three-layer merge (see
[../cli/usage.md](../cli/usage.md#precedence-settingsjson--yaml--flags)).

---

## Top-level identity keys

| Key | Type | Default | Description |
|---|---|---|---|
| `app-name` | string | `"MSpecViz"` | Shown in the Streamlit page title (`st.set_page_config`). |
| `version` | string | `"0.1.0"` | Cosmetic; shown in the sidebar footer. |
| `repository-name` | string | `"MSpecViz"` | Used in default export project names. |
| `online_deployment` | bool | `false` | When `true`, locks down workspace management to prevent path injection; disable for shared-tenant deployments. |
| `enable_workspaces` | bool | `true` | If `false`, the workspace switcher in the sidebar is hidden and all state lives in a single default workspace. |
| `workspaces_dir` | string | `"workspaces"` | Directory (relative to project root) where workspace folders are created. |
| `test` | bool | `false` | Internal flag — set to `true` by the test suite to skip long-running initialisation. |
| `image-format` | string | `"png"` | Image format for Streamlit download buttons that render inline previews. |

---

## `data_directories`

| Key | Type | Default | Description |
|---|---|---|---|
| `data_directories.default` | string | `"/mnt/volume"` | Starting directory for the Setup tab's file browser. Point this at your lab's shared storage root. Symlinks are followed. |

---

## `models`

Catalogue of Koina model IDs the UI exposes. Only models listed here
appear in the Setup tab's Koina dropdown.

| Key | Type | Description |
|---|---|---|
| `models.default` | string | Model pre-selected on first launch. |
| `models.available` | list[string] | Every Koina model the UI will offer. Valid IDs are the ones accepted by the Koina API — see the [Koina model zoo](https://koina.wilhelmlab.org/docs). |

**Shipped defaults** (14 models): `ms2pip_timsTOF2024`, `ms2pip_HCD2021`,
`ms2pip_Immuno_HCD`, `ms2pip_TTOF5600`, `ms2pip_CID_TMT`,
`ms2pip_iTRAQphospho`, `Prosit_2019_intensity`,
`Prosit_2020_intensity_HCD`, `Prosit_2020_intensity_CID`,
`Prosit_2020_intensity_TMT`, `Prosit_2023_intensity_timsTOF`,
`AlphaPept_ms2_generic`, `UniSpec`, `PredFull`.

---

## `koina`

Defaults for AI-prediction settings. Override per-workspace through
`default-parameters.json` and per-run via the Setup tab or YAML.

| Key | Type | Default | Description |
|---|---|---|---|
| `koina.enabled` | bool | `true` | Master switch — if `false`, no Koina fallback is attempted even when a query has no reference. |
| `koina.model` | string | `"ms2pip_timsTOF2024"` | Model used when a query has no experimental reference. Must be in `models.available`. |
| `koina.collision_energy` | int | `35` | Normalized collision energy (NCE) passed to the Koina API. Typical range: 20–40 for HCD/CID. |
| `koina.instrument_type` | string | `"TIMSTOF"` | Instrument hint. Allowed: `QE`, `LUMOS`, `TIMSTOF`, `SCIEXTOF`. Model-dependent — UniSpec and PredFull ignore it. |

---

## `ion_types`

| Key | Type | Default | Description |
|---|---|---|---|
| `ion_types.default` | list[string] | `["b", "y"]` | Fragment ion series pre-selected for new workspaces. |
| `ion_types.available` | list[string] | `["a", "b", "c", "x", "y", "z", "precursor", "immonium"]` | Full set of ion series the UI offers. Values outside this list are ignored by the annotation step. |

---

## `mz_tolerance`

| Key | Type | Default | Description |
|---|---|---|---|
| `mz_tolerance.default` | float | `0.02` | Default m/z matching tolerance in Da. |
| `mz_tolerance.unit` | string | `"Da"` | Informational — the annotation code always treats `mz_tolerance` as Daltons. |

---

## `annotation`

| Key | Type | Default | Description |
|---|---|---|---|
| `annotation.add_losses` | bool | `false` | Include neutral-loss ions (H₂O and NH₃) in the theoretical fragment list. Roughly doubles the number of expected m/z values. |
| `annotation.add_first_prefix_ion` | bool | `true` | Include a1/b1/c1 (the first prefix ion of each series). Useful for short peptides; otherwise the series starts at index 2. |
| `annotation.remove_precursor_peaks` | bool | `false` | Strip peaks within `mz_tolerance` of the precursor m/z before annotating. Reduces noise from unfragmented precursor. |

---

## `performance`

Global tuning. Defaults are conservative; increase `workers` when you run
on a box with spare cores.

| Key | Type | Default | Description |
|---|---|---|---|
| `performance.cache_size_mb` | int | `500` | Parquet cache budget per workspace. Oldest entries are evicted first. |
| `performance.max_file_size_gb` | int | `3` | Upper limit on any single mzML. Files larger than this are skipped with a warning to avoid OOM. |
| `performance.chunk_size` | int | `5000` | Number of spectra streamed before yielding back to the event loop in Streamlit. Increase for throughput on large mzMLs, decrease to keep progress bars snappy. |
| `performance.workers` | int | `4` | Default parallelism for `search_peptides_in_idxmls`, `load_spectra_streaming`, and `run_parallel_comparisons`. `1` = sequential fallback (useful for debugging). |
| `performance.tmpfs_staging` | bool | `true` | Copy mzML files into a tmpfs path before random-access indexing — massively speeds up random access when the source FS is slow/networked. |
| `performance.tmpfs_path` | string | `"/dev/shm"` | tmpfs mount point. |
| `performance.tmpfs_ram_fraction` | float | `0.40` | Maximum fraction of system RAM to consume for tmpfs staging. Files above this cap fall back to direct FS access. |

---

## `export`

Global export defaults; overridable per-workspace and per-CLI-run.

| Key | Type | Default | Description |
|---|---|---|---|
| `export.default_directory` | string | `"./exports"` | Base directory for versioned export subdirs (see [../export/layout.md](../export/layout.md)). |
| `export.project_name` | string | `"MSpecViz_analysis"` | Default project name prefix (e.g. `MSpecViz_analysis_001_20260415/`). |
| `export.formats` | list[string] | `["png"]` | Default set of static image formats. Allowed: `png`, `svg`, `pdf`. HTML (interactive Plotly) is a separate export flag. |
| `export.dpi` | int | `300` | Raster DPI for PNG outputs. |
| `export.figsize` | list[int, int] | `[12, 8]` | Figure size in inches `[width, height]`. |

---

## `cache`

| Key | Type | Default | Description |
|---|---|---|---|
| `cache.enabled` | bool | `true` | When `true`, loaded spectra + annotated DataFrames are persisted to `<workspace>/.mspecviz_cache/` keyed by settings fingerprint. Re-runs with identical settings skip the load stage. |
| `cache.subdirectory` | string | `".mspecviz_cache"` | Directory name under each workspace. Leading dot makes it hidden by default in file browsers. |

---

## `ion_colors`

Colour palette used by both mirror renderers. Accepts any string Plotly +
Matplotlib both understand (hex, named CSS colours, `rgb()` strings).

| Key | Default (hex) | Used for |
|---|---|---|
| `ion_colors.a` | `#FFCC00` | a-ion series |
| `ion_colors.b` | `#0066CC` | b-ion series |
| `ion_colors.c` | `#00CC66` | c-ion series |
| `ion_colors.x` | `#CC6600` | x-ion series |
| `ion_colors.y` | `#CC0000` | y-ion series |
| `ion_colors.z` | `#6600CC` | z-ion series |
| `ion_colors.precursor` | `#8B008B` | Precursor ion and its charge states |
| `ion_colors.immonium` | `#00CED1` | Immonium ions |
| `ion_colors.unannotated` | `#AAAAAA` | Peaks with no theoretical match |

Charge state (`++`, `+++`) does not change the base colour. The palette
is exported as `ION_COLORS` from `mspecviz.plotting` and also used in
summary plots.
