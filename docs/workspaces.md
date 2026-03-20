# Workspaces

A **workspace** is a named folder where MSpecViz saves everything from one analysis project — your file selections, settings, result images, and summary tables. Each workspace is completely independent, so different projects never interfere with each other.

---

## What gets saved automatically

Every time you finish a search or export, MSpecViz saves your current settings to the workspace. The next time you open it, everything is exactly as you left it:

- Which idXML and mzML files you selected
- Your peptide list
- Ion types, m/z tolerance, and annotation options
- Export format and project name

You never need to re-enter these manually.

---

## Creating and switching workspaces

Open the **Workspaces** panel in the left sidebar.

- **To create a new workspace:** type a name and click **Create Workspace**. Use letters, numbers, and hyphens — no spaces or slashes.
- **To switch:** choose an existing workspace from the dropdown. Your previous settings load automatically.
- **To delete:** type the name and click **Delete Workspace**. This permanently removes all files inside it. The `default` workspace cannot be deleted.

---

## What's inside a workspace folder

After running an analysis, your workspace folder on disk looks like this:

```
my_experiment/
├── params.json                    ← your saved settings
└── exports/
    └── MSpecViz_analysis_20260313/
        ├── png/                   ← mirror plot images
        ├── svg/
        ├── pdf/
        └── MSpecViz_analysis_summary.xlsx
```

The `exports/` folder contains all your plots and the summary Excel file. Each analysis run gets its own dated sub-folder, so older results are never overwritten.

---

## Sharing results with a colleague

Because a workspace is just a folder, sharing is simple:

**Option 1 — Share only the results** (plots + Excel summary):
Copy the `exports/` sub-folder and send it by email, shared drive, or any file transfer method. No MSpecViz installation required on the recipient's side.

**Option 2 — Share the full workspace** (results + settings):
Copy the entire workspace folder to the recipient's `workspaces-MSpecViz/` directory. When they open MSpecViz, the workspace will appear in their sidebar dropdown. Selecting it restores your exact settings and shows all previous export runs in the History tab.

> **Tip:** The hidden `.mspecviz_cache/` sub-folder can be large (it speeds up re-runs by caching the raw data). You can safely skip it when sharing — MSpecViz will rebuild it automatically if needed.

---

## Browsing previous runs

The **8. History** tab shows all past export runs for the current workspace. For each run you can:

- See the date, peptide list, and file counts at a glance
- Click **Reload settings** to restore that run's settings into the current session
- Click **Delete run** to remove the export folder from disk (with confirmation)

This is useful when you want to repeat an analysis with slightly different settings, or when a colleague has shared their workspace and you want to start from their configuration.
