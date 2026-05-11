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
YAML and editing it. The log level is controlled by an environment
variable:

```bash
# Default: INFO (stage headers, per-file logs, progress bars)
python -m cli my_config.yaml

# Quiet — errors only
MSPECVIZ_LOG_LEVEL=ERROR python -m cli my_config.yaml

# Verbose — full DEBUG from pyopenms, koinapy, mspecviz.*
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

1. **`mspecviz/common/defaults.py`** seeds the runtime config with every
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
| `1` | Pipeline failed (pyopenms error, Koina unreachable, filesystem I/O, …). Error is logged and the process exits. |
| `2` | Config error (YAML path missing, YAML not a mapping, required section missing, peptide list empty). |

The test suite uses `test_missing_config_exits_two` to lock the
exit-code behaviour, so it's safe to automate against.

---

## Progress output

Each stage logs a header line, and a `tqdm` progress bar tracks
file/spectrum/plot completion within the stage. Set
`MSPECVIZ_LOG_LEVEL=ERROR` to suppress everything except real errors;
the progress bars stay on by default (they go to stderr).

Typical successful run:

```text
[INFO] Step 1: Searching peptides in idXML files
Searching: 100%|██████████| 1/1 [00:01<00:00,  1.2file/s]
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
