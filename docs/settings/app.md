# App-wide settings — `settings.json`

Path: `streamlit_app/config/settings.json`

App-wide environmental configuration and enumerations. Read once at
startup into `st.session_state.settings`; every value is looked up via
`get_nested(settings, "section", "subkey", default=...)` so a missing
key falls back to a hard-coded default.

`settings.json` does **not** own tunable-parameter defaults. Those live
in `core/common/defaults.py` (`DEFAULTS`) — see
[overview.md](overview.md) and
[../../core/common/defaults.py](../../core/common/defaults.py).

This file is Streamlit-only. The CLI reads `DEFAULTS` directly and does
not open `settings.json`.

---

## Top-level identity keys

| Key | Type | Default | Description |
|---|---|---|---|
| `app_name` | string | `"MSpecViz"` | Shown in the Streamlit page title (`st.set_page_config`). |
| `version` | string | `"0.1.0"` | Cosmetic; shown in the sidebar footer. |
| `repository_name` | string | `"MSpecViz"` | Used in default export project names. |
| `online_deployment` | bool | `false` | When `true`, locks down workspace management to prevent path injection; disable for shared-tenant deployments. |
| `enable_workspaces` | bool | `true` | If `false`, the workspace switcher in the sidebar is hidden and all state lives in a single default workspace. |
| `workspaces_dir` | string | `"workspaces"` | Directory (relative to project root) where workspace folders are created. |
| `test` | bool | `false` | Internal flag — set to `true` by the test suite to skip long-running initialisation. |

---

## `data_directories`

| Key | Type | Default | Description |
|---|---|---|---|
| `data_directories.default` | string | `"/mnt/volume"` | Starting directory for the Setup tab's file browser. Point this at your lab's shared storage root. Symlinks are followed. |

---

## `models`

Catalogue of Koina model IDs the UI exposes. Only models listed here
appear in the Setup tab's Koina dropdown. The pre-selected model comes
from `DEFAULTS.koina_model`, not this file.

| Key | Type | Description |
|---|---|---|
| `models.available` | list[string] | Every Koina model the UI will offer. Valid IDs are the ones accepted by the Koina API — see the [Koina model zoo](https://koina.wilhelmlab.org/docs). |

**Shipped list** (14 models): `ms2pip_timsTOF2024`, `ms2pip_HCD2021`,
`ms2pip_Immuno_HCD`, `ms2pip_TTOF5600`, `ms2pip_CID_TMT`,
`ms2pip_iTRAQphospho`, `Prosit_2019_intensity`,
`Prosit_2020_intensity_HCD`, `Prosit_2020_intensity_CID`,
`Prosit_2020_intensity_TMT`, `Prosit_2023_intensity_timsTOF`,
`AlphaPept_ms2_generic`, `UniSpec`, `PredFull`.

---

## `ion_types`

| Key | Type | Description |
|---|---|---|
| `ion_types.available` | list[string] | Full set of ion series the UI offers. Values outside this list are ignored by the annotation step. |

Shipped: `["a", "b", "c", "x", "y", "z", "internal", "precursor",
"immonium"]`. The pre-selected ion types come from `DEFAULTS.ion_types`
(which defaults to `("b", "y")`). `internal` corresponds to
spectrum_utils' "m" ions (internal fragments).

---

## `mz_tolerance`

| Key | Type | Default | Description |
|---|---|---|---|
| `mz_tolerance.unit` | string | `"Da"` | Informational — the annotation code always treats `mz_tolerance` as Daltons. |

Default numeric value lives in `DEFAULTS.mz_tolerance`.

---

## `performance`

Global tuning for environmental knobs (cache sizing, tmpfs staging).
Worker count is a pipeline parameter and lives in `DEFAULTS.workers`.

| Key | Type | Default | Description |
|---|---|---|---|
| `performance.cache_size_mb` | int | `500` | Parquet cache budget per workspace. Oldest entries are evicted first. |
| `performance.max_file_size_gb` | int | `3` | Upper limit on any single mzML. Files larger than this are skipped with a warning to avoid OOM. |
| `performance.chunk_size` | int | `5000` | Number of spectra streamed before yielding back to the event loop in Streamlit. Increase for throughput on large mzMLs, decrease to keep progress bars snappy. |
| `performance.tmpfs_staging` | bool | `true` | Copy mzML files into a tmpfs path before random-access indexing — massively speeds up random access when the source FS is slow/networked. |
| `performance.tmpfs_path` | string | `"/dev/shm"` | tmpfs mount point. |
| `performance.tmpfs_ram_fraction` | float | `0.40` | Maximum fraction of system RAM to consume for tmpfs staging. Files above this cap fall back to direct FS access. |

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

The shipped palette is Material Design 700-level (saturated, print-safe):

| Key | Default (hex) | Used for |
|---|---|---|
| `ion_colors.a` | `#388E3C` | a-ion series (Material Green 700) |
| `ion_colors.b` | `#1976D2` | b-ion series (Material Blue 700) |
| `ion_colors.c` | `#00796B` | c-ion series (Material Teal 700) |
| `ion_colors.x` | `#7B1FA2` | x-ion series (Material Purple 700) |
| `ion_colors.y` | `#D32F2F` | y-ion series (Material Red 700) |
| `ion_colors.z` | `#F57C00` | z-ion series (Material Orange 700) |
| `ion_colors.internal` | `#FBC02D` | Internal fragments / spectrum_utils' "m" ions (Material Yellow 700) |
| `ion_colors.precursor` | `#512DA8` | Precursor ion and its charge states (Material Deep Purple 700) |
| `ion_colors.immonium` | `#455A64` | Immonium ions (Material Blue Grey 700) |
| `ion_colors.unannotated` | `#9E9E9E` | Peaks with no theoretical match (Material Grey 500) |

Charge state (`++`, `+++`) does not change the base colour. The same
palette is the single source of truth: it ships in `settings.json` for
the Streamlit UI and as `ION_COLORS` in
`core/common/defaults.py` for everything else (mirror renderer,
summary plots, CLI). Editing one without the other will desynchronise
the two — keep them in lockstep or remove the JSON copy.
