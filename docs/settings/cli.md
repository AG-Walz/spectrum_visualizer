# CLI Settings — YAML & Flags

The two configuration surfaces specific to `python -m cli`:

- a YAML config file (one per run, passed via `-c`)
- argparse command-line flags (single-invocation overrides)

Every YAML key and every CLI flag `dest=` matches a field name on
`CLIConfig` exactly — no translation table. Defaults come from
`mspecviz.common.defaults.DEFAULTS`. See [overview.md](overview.md) for
precedence; [../cli/usage.md](../cli/usage.md) for recipes.

---

## YAML config structure

Flat snake_case keys — one name per concept, matching `CLIConfig`
fields one-for-one. Every key is optional except `idxml_files`,
`mzml_files`, and `peptides`.

```yaml
# Inputs
idxml_files: []                       # list of paths or glob patterns
mzml_files:  []                       # list of paths or glob patterns
peptide_source: all                   # "all" | /path/to/peptides.txt | [PEP1, PEP2]
reference_substring: ""               # empty = single-spectrum mode

# Annotation
ion_types: [b, y]                     # any of: a b c x y z precursor immonium
mz_tolerance: 0.02                    # Da
add_neutral_losses: false             # H2O + NH3 losses
add_first_prefix_ion: true            # include a1/b1/c1
remove_precursor_peaks: false         # strip precursor before annotating
intensity_mode: Relative              # "Relative" or "Raw"

# Export
export_directory: ./exports           # base output dir
project_name: MSpecViz_analysis       # subdir prefix
export_formats: [png]                 # png svg pdf html (any combo)
export_resolution: 300                # raster DPI (dots per inch)
figure_size: [12, 8]                  # inches [width, height]
plot_unpaired_queries: false          # emit single-spectrum plots when no ref

# Performance
workers: 4                            # 1 = sequential fallback

# Summary plots
summary_metric: "R² binned"           # default summary-plot metric
summary_plot_types: [dual_axis, ridgeline]  # subset of dual_axis, ridgeline, heatmap

# Koina AI prediction
koina_enabled: false                  # master switch for AI prediction
koina_model: ms2pip_timsTOF2024       # any model in settings.json > models.available
koina_collision_energy: 35            # NCE
koina_instrument_type: TIMSTOF        # QE | LUMOS | TIMSTOF | SCIEXTOF
```

---

## Key reference

| Key | Default | Description |
|---|---|---|
| `idxml_files` | `[]` | idXML paths or glob patterns. Globs are expanded with `recursive=True`. |
| `mzml_files` | `[]` | mzML paths or glob patterns. |
| `peptide_source` | *(required)* | `"all"`, `/abs/path/peptides.txt`, or inline list `[PEP1, PEP2]`. |
| `reference_substring` | `""` | Filename substring identifying reference mzMLs. Empty ⇒ single-spectrum mode. |
| `ion_types` | `[b, y]` | Fragment ion series. |
| `mz_tolerance` | `0.02` | Da. |
| `add_neutral_losses` | `false` | Neutral losses (H₂O, NH₃). |
| `add_first_prefix_ion` | `true` | Include a1/b1/c1. |
| `remove_precursor_peaks` | `false` | Strip precursor peaks before annotation. |
| `intensity_mode` | `"Relative"` | `"Relative"` or `"Raw"`. |
| `export_directory` | `"./exports"` | Base export dir. |
| `project_name` | `"MSpecViz_analysis"` | Directory prefix. |
| `export_formats` | `[png]` | Any of `png`, `svg`, `pdf`, `html`. Multiple formats can be combined. |
| `export_resolution` | `300` | Raster DPI (dots per inch) for PNG/PDF. |
| `figure_size` | `[12, 8]` | Inches [width, height]. |
| `plot_unpaired_queries` | `false` | When `true`, single-spectrum plots are generated for queries with no reference. When `false`, those queries are skipped. |
| `workers` | `4` | Parallel workers. `1` forces the sequential fallback (useful for debugging). |
| `summary_metric` | `"R² binned"` | Metric used to colour/sort summary plots. Must be one of the keys in `METRIC_COLUMNS`. |
| `summary_plot_types` | `[dual_axis, ridgeline]` | Subset of `{dual_axis, ridgeline, heatmap}` to render. |
| `koina_enabled` | `false` | Master switch. `--no-koina` flips it off. |
| `koina_model` | `"ms2pip_timsTOF2024"` | Model ID. |
| `koina_collision_energy` | `35` | NCE. |
| `koina_instrument_type` | `"TIMSTOF"` | Instrument hint. |

### `peptides` interpretation

| Value | Behaviour |
|---|---|
| `"all"` | Scan every selected idXML and use every unique peptide (capped at 10 000). |
| `/abs/path/peptides.txt` | Read newline- or comma-separated peptide list from disk. |
| `[PEPTIDE1, PEPTIDE2]` | Inline list (YAML flow style). |

Peptides are case-folded to uppercase and stripped of whitespace before
the search.

---

## Command-line flags

Every flag overrides the corresponding YAML key for a single invocation.
The argparse `dest=` matches the `CLIConfig` field name, so the mapping
is direct.

| Flag | Argument | Overrides |
|---|---|---|
| `-c`, `--config` | PATH | *(required)* path to YAML config |
| `--idxml` | PATTERN… | `idxml_files` |
| `--mzml` | PATTERN… | `mzml_files` |
| `--peptides` | PATH_OR_ALL | `peptides` |
| `--reference-substring` | STR | `reference_substring` |
| `--ion-types` | ION… | `ion_types` |
| `--mz-tolerance` | DA | `mz_tolerance` |
| `--export-directory` | DIR | `export_directory` |
| `--project-name` | NAME | `project_name` |
| `--formats` | FMT… | `export_formats` (choices: png svg pdf html) |
| `--workers` | N | `workers` |
| `--summary-metric` | NAME | `summary_metric` |
| `--summary-plot-types` | TYPE… | `summary_plot_types` (choices: dual_axis ridgeline heatmap) |
| `--no-plots` | *(flag)* | Skip plot generation; export CSV/Parquet only |
| `--no-koina` | *(flag)* | Force `koina_enabled = false` regardless of YAML |
| `--koina-model` | ID | `koina_model` |
| `--koina-collision-energy` | EV | `koina_collision_energy` |
| `--koina-instrument-type` | TYPE | `koina_instrument_type` |
| `-v`, `--verbose` | *(flag)* | Enable DEBUG logging |
| `-q`, `--quiet` | *(flag)* | Suppress non-error output |

---

## See also

- [../cli/usage.md](../cli/usage.md) — invocation, exit codes, recipes.
- [app.md](app.md) — the `settings.json` layer (Streamlit-only).
- [overview.md](overview.md) — precedence rules.
