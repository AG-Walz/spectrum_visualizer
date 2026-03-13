# Workspaces

MSpecViz uses a **workspace** system to keep each project's files, exports, and cache
isolated from one another.  This document explains how workspaces work, how to manage
them, and how to share them with colleagues.

---

## What is a workspace?

A workspace is a named directory on disk that stores everything produced by one
analysis session:

| Sub-directory | Contents |
|---------------|----------|
| `exports/<project>_<hash>/` | Versioned export folders with mirror-plot images (PNG/SVG/PDF) and the summary Excel/CSV file |
| `.mspecviz_cache/` | Parquet cache for fast re-runs without re-loading mzML/idXML files |
| `params.json` | Last-used settings for that workspace (restored on next visit) |

Workspaces live inside the parent directory configured by `workspaces_dir` in
`settings.json` (default: `..`) under a project-specific sub-directory:

```
<workspaces_dir>/
└── workspaces-MSpecViz/
    ├── default/         ← the built-in workspace, always present
    ├── my_experiment/
    └── validation_run/
```

By default, workspaces are stored one level above the application directory:

```
parent/
├── MSpecViz/            ← application code
└── workspaces-MSpecViz/ ← all workspace directories
```

You can change `workspaces_dir` in `settings.json` to any absolute or relative path.

---

## Enabling workspaces

Workspaces are enabled by default.  The relevant keys in `settings.json` are:

```json
{
  "enable_workspaces": true,
  "workspaces_dir": ".."
}
```

Set `"enable_workspaces": false` to disable the switcher UI entirely; all data will
still be written to the `default` workspace directory.

---

## Creating a workspace (UI)

1. Open the **Workspaces** panel in the left sidebar (the expander labelled
   **Workspaces**).
2. Type a new name (letters, numbers, hyphens, and underscores — no slashes) into the
   **create/remove workspace** field.
3. Click **Create Workspace**.

The app switches to the new workspace immediately.  Its directory is created on disk
and will appear in the workspace selector dropdown on your next visit.

---

## Switching between workspaces (UI)

Use the **choose existing workspace** dropdown in the sidebar **Workspaces** panel.
Selecting a workspace reloads the page and restores its `params.json` settings.

The active workspace is also encoded in the browser URL as a query parameter, e.g.:

```
http://localhost:8501/?workspace=my_experiment
```

You can bookmark this URL to jump directly to a specific workspace.

---

## Deleting a workspace (UI)

1. Type the workspace name into the **create/remove workspace** field.
2. Click **Delete Workspace**.

The entire workspace directory (including all exports and cache) is removed from disk.
The `default` workspace cannot be deleted.

---

## Directory structure after a run

After processing, a workspace looks like this:

```
workspaces-MSpecViz/my_experiment/
├── params.json                          ← saved settings
├── exports/
│   └── MSpecViz_analysis_20260313_141500_a1b2c3d4/
│       ├── png/                         ← mirror plot images
│       │   ├── PEPTIDE_z2_run1.png
│       │   └── ...
│       ├── svg/
│       ├── pdf/
│       └── MSpecViz_analysis_summary.xlsx
└── .mspecviz_cache/
    └── <fingerprint>.parquet            ← re-run cache
```

Each export subfolder is **fingerprinted** (SHA-256 of all run settings) so plots from
different parameter combinations never overwrite each other.

---

## Sharing a workspace

Because a workspace is just a directory, sharing it is straightforward:

1. **Archive** the workspace directory:
   ```bash
   tar -czf my_experiment.tar.gz workspaces-MSpecViz/my_experiment/
   ```
2. **Transfer** the archive (SCP, shared network drive, object storage, etc.).
3. The recipient places it in their `workspaces-MSpecViz/` directory and selects it
   from the sidebar dropdown.

> **Tip:** The Parquet cache in `.mspecviz_cache/` can be large.  Exclude it if you
> only want to share the final plots and summary:
> ```bash
> tar -czf my_experiment_results.tar.gz \
>   --exclude='my_experiment/.mspecviz_cache' \
>   workspaces-MSpecViz/my_experiment/
> ```

---

## Workspace in the CLI

The CLI (`cli.py`) does not interact with workspaces directly — you specify the export
directory explicitly:

```yaml
# config.yaml
export:
  directory: /path/to/workspaces-MSpecViz/my_experiment/exports
```

or as a CLI flag:

```bash
python cli.py --config config.yaml --export-dir /path/to/workspaces-MSpecViz/my_experiment/exports
```

Use the same `exports/` sub-directory as the Streamlit app to keep CLI and GUI outputs
in the same workspace.

---

## Template background

The workspace design follows the
[OpenMS streamlit-template](https://github.com/OpenMS/streamlit-template) convention.
`page_setup()` in `src/common/common.py` initialises `st.session_state.workspace` at
application startup and is called at the top of every page.  The workspace path is then
available throughout the app as `st.session_state.workspace`.
