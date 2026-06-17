# Settings

Where configuration comes from and how the layers merge. The field list with
default values is **not** here — it lives once in `core/defaults.py` (`DEFAULTS`,
a frozen dataclass). Every other layer only overrides it.

---

## Sources and precedence

| Source | Scope | Read by |
|---|---|---|
| `core/defaults.py` (`DEFAULTS`) | Tunable-parameter defaults (single source of truth) | Both front ends |
| `streamlit_app/config/settings.json` | App-wide env config + enumerations | Streamlit only |
| `workspaces/<name>/params.json` | One workspace's edited widget values | Streamlit only |
| `cli/config/*.yaml` | One CLI invocation | CLI only (see [cli.md → YAML config](cli.md#yaml-config)) |
| `MSPECVIZ_LOG_LEVEL` env var | One invocation | Both — logging verbosity |

All paths resolve against the **project root** (the Streamlit app `cd`s into
`streamlit_app/` only so it picks up `.streamlit/config.toml`).

- **CLI run:** `DEFAULTS < YAML`. `cli/main.py:_load_config` seeds from `DEFAULTS`
  (plus CLI-only inputs), overlays the YAML by exact field name, no flag overrides.
- **Streamlit run:** `DEFAULTS < workspace params.json < live widget state`.
  `settings.json` is read in parallel for enumerations and env flags — it does
  **not** carry tunable defaults.

**Naming convention.** One `snake_case` name per setting, identical across every
layer (the `DEFAULTS` field, the YAML key, the `params.json` key, the session-state
key). No aliases, no kebab/camelCase. A missing key always falls back to
`DEFAULTS.<field>` at the call site. To change a default, edit `core/defaults.py`;
that affects every workspace and CLI run that hasn't overridden it.

---

## `settings.json` (Streamlit-only)

App-wide environmental config, read once into `st.session_state.settings`. It holds
only env-level config the defaults can't express — the app name, the workspace
toggle and location, the file-browser start directory, and the full ion-series set
the UI offers (the pre-selected subset comes from `DEFAULTS.ion_types`). It does
**not** own tunable defaults, and the Koina catalogue (`KOINA_MODELS`) and ion
palette (`ION_COLORS`) live in `core/defaults.py`, not here. The file is short and
self-documenting — read `streamlit_app/config/settings.json` for the current keys.

---

## Workspaces and session state

A new workspace starts with an empty `params.json` — there is no template file.
`utils/state.py:save_params()` writes a key only when the user changes its widget;
every read falls back to `DEFAULTS.<key>` (via `get_param`). Names are validated by
`is_safe_workspace_name()` (no `..`, absolute paths, or shell metacharacters).

Most `params.json` keys are tunables already defined in `DEFAULTS`. The exceptions
are the **Setup-tab inputs**, which have no other home: the queued idXML/mzML paths,
the last-browsed directory per picker, the raw peptide textarea, and the reference
substring (empty ⇒ single-spectrum mode).

Every widget writes into `st.session_state` under the `peptide_mirror.` namespace
(`STATE_PREFIX`, a UI-layer concern — always resolve via `get_state_key(key)`).
Beyond those parameters, the keys that matter as a **cross-view contract** are the
derived artifacts Load & Annotate produces and Explore / Export consume: the four
Polars LazyFrames, the PSM lookup, the annotated queries, the reference store, the
comparison frame, and the `peptidoform_psms` grouping — gated by a `data_ready` flag
and a `cache_version` counter bumped on each fresh load. The authoritative list is
the set persisted to the analysis cache in `utils/page_state.py` (`_LAZYFRAME_KEYS`,
`_PICKLE_KEYS`, and the explicit handling around them); because every view reads its
state back through that cache, the contract and its serialiser change together.
