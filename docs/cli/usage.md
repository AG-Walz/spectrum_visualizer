# CLI Usage

`python -m cli` runs the full MSpecViz pipeline headlessly, reading its
inputs from a YAML config and writing results to a versioned export
directory — same output format as the Streamlit Export tab.

For the full schema of every flag and YAML key, see
[../settings/cli.md](../settings/cli.md). For recipes and debugging tips,
see [recipes.md](recipes.md).

---

## Invocation

From the project root:

```bash
python -m cli -c cli/config/config_example.yaml
```

The `-c` / `--config` flag is the only required argument. Everything else
has a default in either the YAML file or `streamlit_app/config/settings.json`.

Useful one-liners:

```bash
# Verbose (DEBUG-level logs)
python -m cli -c my_config.yaml -v

# Quiet (errors only)
python -m cli -c my_config.yaml -q

# Override workers for a single run
python -m cli -c my_config.yaml --workers 8

# Export only DataFrames + CSV (no plots)
python -m cli -c my_config.yaml --no-plots

# Force Koina off regardless of config
python -m cli -c my_config.yaml --no-koina

# Custom output directory
python -m cli -c my_config.yaml --export-directory /tmp/my_exports
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

## Precedence: `DEFAULTS` < YAML < flags

Two override layers, merged by `cli/config_builder.py`:

1. **`mspecviz/common/defaults.py`** seeds `CLIConfig` with every
   default.
2. **YAML file** overrides any seeded value.
3. **Command-line flags** override both above.

YAML keys and argparse `dest=` names match `CLIConfig` field names
exactly, so the merge is a direct name-for-name copy — no translation
table. To change a default for all CLI runs, edit `defaults.py`. To
change it per run, use YAML. To change it for one invocation, use a
flag. Adding a new parameter is a one-field edit on `CLIConfig`.

---

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Pipeline completed successfully. |
| `1` | Pipeline failed (pyopenms error, Koina unreachable, filesystem I/O, …). Error is logged and the process dies. |
| `2` | Config error (YAML path missing, YAML not a mapping, required section missing, unknown format in `--formats`). |

The test suite uses `test_missing_config_exits_two` to lock down the
exit-code-2 behaviour, so it's safe to automate against.

---

## Progress output

By default the CLI prints a stage header for each pipeline step and a
terminal progress bar (`cli.output.TerminalProgress`). With `-q`, only
stage headers and final error lines are printed. With `-v`, DEBUG logs
from `pyopenms`, `koinapy`, and internal helpers appear inline.

Typical successful run:

```
Stage 1/3: Search peptides
  Searching 1 idXML file... done (2 variants)
Stage 2/3: Load spectra
  Streaming 8 mzML files... done (2 queries, 0 refs, 2 Koina fallbacks)
Stage 3/3: Build DataFrames + export
  PSM table: 2 rows
  Peaks: 373 rows
  Writing run_settings.json... done
  Exporting 2 mirror plots (png)... done

Pipeline Complete:
  Variants: 2 | PSMs: 2 | Spectra: 2 | Comparisons: 2 | Elapsed: 17s
  Output: exports/MSpecViz_analysis_001_20260415/
```

---

## Output layout

The CLI writes to the **same versioned directory structure** as the
Streamlit Export tab — see [../export/layout.md](../export/layout.md).
Briefly:

```
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
