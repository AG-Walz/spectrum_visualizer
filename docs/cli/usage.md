# CLI Usage

`python -m cli` runs the full MSpecViz pipeline headlessly, reading its
inputs from a YAML config and writing results to a versioned export
directory — same output format as the Streamlit Export tab.

For the YAML schema, see [../settings/cli.md](../settings/cli.md). For
recipes and debugging tips, see [recipes.md](recipes.md).

---

## Invocation

From the project root:

```bash
python -m cli cli/config/config_example.yaml
```

A single positional argument: the path to a YAML config file. There are
no command-line flag overrides — per-run tweaks are made by copying the
YAML and editing it. Path globs in `idxml_files` / `mzml_files` are
expanded by `glob.glob(pattern, recursive=True)`, so `**/` wildcards
work (`/data/**/run_*.mzML`).

The log level is controlled by an environment variable:

```bash
# Default: INFO (stage headers, per-file logs, progress bars)
python -m cli my_config.yaml

# Quiet — errors only
MSPECVIZ_LOG_LEVEL=ERROR python -m cli my_config.yaml

# Verbose — full DEBUG from pyopenms, koinapy, core.*
MSPECVIZ_LOG_LEVEL=DEBUG python -m cli my_config.yaml
```

---

## Minimum viable config

Smallest YAML that will run (everything else inherits from `DEFAULTS`):

```yaml
idxml_files:
  - /path/to/results.idXML
mzml_files:
  - /path/to/raw/*.mzML

peptide_source: /path/to/peptides.txt
```

The sample at `cli/config/config_example.yaml` is a fully-annotated
template — copy it and edit the paths.

---

## Precedence: `DEFAULTS` < YAML

Two layers, merged inside `cli/main.py:_load_config`:

1. **`core/defaults.py`** seeds the runtime config with every
   default plus CLI-only fields (`idxml_files`, `mzml_files`,
   `peptide_source`, `reference_substring`, `no_plots`).
2. **YAML file** overrides any seeded value.

YAML keys match the field names on the merged config one-for-one — no
translation table. To change a default for all CLI runs, edit
`defaults.py`. To change it per run, edit the YAML.

---

## Exit codes

| Code | Meaning |
| --- | --- |
| `0` | Pipeline completed successfully. |
| Non-zero | Anything went wrong. |

The CLI uses `sys.exit("Error: …")` for every failure path (missing
config, malformed YAML, missing inputs, no PSMs, no spectra, or a
pyopenms/Koina/filesystem exception in the pipeline). That always yields
exit code **1** with the message on stderr — there is no separate
"config error" code. Test `test_missing_config_exits_two` (in
`tests/cli/test_cli_integration.py`) only asserts non-zero despite its
name; do not rely on the exact code for branching.

---

## Progress output

Each stage logs a header line, and a `tqdm` progress bar tracks
file/spectrum/plot completion within the stage. The bars are driven by
the same `progress_callback` / `error_callback` hooks the Streamlit
front end uses — both layers feed updates into
`core.helpers.run_tasks` and let each frontend render
them its own way (tqdm vs `st.progress`).

Set `MSPECVIZ_LOG_LEVEL=ERROR` to suppress everything except real
errors; the progress bars stay on by default (they go to stderr).

`workers` (default 4) controls parallelism for all three pipeline
stages — idXML lookup, mzML streaming, and plot rendering. `workers: 1`
forces the sequential code path inside `run_tasks`, which is the right
setting when you want a worker-side traceback instead of a generic
"worker crashed" message.

Typical successful run:

```text
[INFO] Step 1: Looking up peptides in idXML files
Looking up: 100%|██████████| 1/1 [00:01<00:00,  1.2file/s]
[INFO] Found 2 PSMs across 2 variants
[INFO] Step 2: Loading and annotating spectra
Loading:   100%|██████████| 2/2 [00:14<00:00,  7.1s/spec]
[INFO] Step 3: Exporting results
Plots:     100%|██████████| 2/2 [00:02<00:00,  0.9plot/s]
[INFO] Done in 17.3s. Exported to exports/MSpecViz_analysis_001_20260415
```

---

## Output layout

The CLI writes to the **same versioned directory structure** as the
Streamlit Export tab — see [../export/layout.md](../export/layout.md).
Briefly:

```text
{export_dir}/{project_name}_{NNN}_{YYYYMMDD}/
    run_settings.json            # canonical snapshot of this run
    comparison_summary.csv       # 1 row per query-reference comparison
    comparison_summary.xlsx
    psm.parquet                  # PSM table
    peaks.parquet                # peaks
    seq.parquet                  # sequences
    heatmap.parquet              # RT-by-mz
    summary/                     # dual-axis + ridgeline plots
    png/ / svg/ / pdf/           # static mirror plots, one subdir per peptide
    html/                        # interactive Plotly mirror plots
```

If a previous run with identical settings exists, the CLI **reuses the
directory** and skips already-generated plots (same behaviour as the
Streamlit Export tab). A bumped `_last_run` timestamp in
`run_settings.json` marks the re-use.
