# CLI Settings — YAML & Flags

The two configuration surfaces specific to `python -m cli`:

- a YAML config file (one per run, passed via `-c`)
- argparse command-line flags (single-invocation overrides)

See [overview.md](overview.md) for how these interact with
`settings.json`. For invocation recipes and debugging tips, see
[../cli/usage.md](../cli/usage.md).

---

## YAML config structure

Full key list with defaults. Every key is optional except `input.*` and
`peptides:`. Defaults come from `streamlit_app/config/settings.json`
unless noted otherwise.

```yaml
input:
  idxml_files: []                  # list of paths or glob patterns
  mzml_files:  []                  # list of paths or glob patterns

peptides: all                      # "all" | /path/to/peptides.txt | [PEP1, PEP2]

reference:
  substring: ""                    # empty = single-spectrum mode

annotation:
  ion_types: [b, y]                # any of: a b c x y z precursor immonium
  mz_tolerance: 0.02               # Da
  add_losses: false                # H2O + NH3 losses
  add_first_prefix_ion: true       # include a1/b1/c1
  remove_precursor_peaks: false    # strip precursor before annotating
  intensity_mode: Relative         # "Relative" or "Raw"

export:
  directory: ./exports             # base output dir
  project_name: MSpecViz_analysis  # subdir prefix
  formats: [png]                   # png svg pdf html (any combo)
  dpi: 300                         # raster DPI
  figsize: [12, 8]                 # inches [width, height]
  plot_unpaired: false             # emit single-spectrum plots when no ref

performance:
  workers: 4                       # 1 = sequential fallback

koina:
  enabled: true                    # master switch for AI prediction
  model: ms2pip_timsTOF2024        # any model in settings.json > models.available
  collision_energy: 35             # NCE
  instrument_type: TIMSTOF         # QE | LUMOS | TIMSTOF | SCIEXTOF
```

---

## Section-by-section reference

### `input:`

| Key | Type | Description |
|---|---|---|
| `input.idxml_files` | list[string] | idXML paths or glob patterns. Globs are expanded with `recursive=True`. |
| `input.mzml_files` | list[string] | mzML paths or glob patterns. |

### `peptides:`

A single value, interpreted one of three ways:

| Value | Behaviour |
|---|---|
| `"all"` | Scan every selected idXML and use every unique peptide (capped at 10 000). |
| `/abs/path/peptides.txt` | Read newline- or comma-separated peptide list from disk. |
| `[PEPTIDE1, PEPTIDE2]` | Inline list (YAML flow style). |

Peptides are case-folded to uppercase and stripped of whitespace before
the search.

### `reference:`

| Key | Type | Description |
|---|---|---|
| `reference.substring` | string | Filename substring identifying reference mzMLs. Empty string ⇒ single-spectrum mode (no mirror plots, no metric comparisons). |

### `annotation:`

Mirrors `settings.json > annotation` plus two extras that are YAML-only:

| Key | Default | Description |
|---|---|---|
| `annotation.ion_types` | `[b, y]` | Fragment ion series. |
| `annotation.mz_tolerance` | `0.02` | Da. |
| `annotation.add_losses` | `false` | Neutral losses (H₂O, NH₃). |
| `annotation.add_first_prefix_ion` | `true` | Include a1/b1/c1. |
| `annotation.remove_precursor_peaks` | `false` | Strip precursor peaks before annotation. |
| `annotation.intensity_mode` | `"Relative"` | `"Relative"` or `"Raw"`. CLI-only key. |

### `export:`

| Key | Default | Description |
|---|---|---|
| `export.directory` | `"./exports"` | Base export dir. |
| `export.project_name` | `"MSpecViz_analysis"` | Directory prefix. |
| `export.formats` | `[png]` | Any of `png`, `svg`, `pdf`, `html`. Multiple formats can be combined. |
| `export.dpi` | `300` | Raster DPI. |
| `export.figsize` | `[12, 8]` | Inches [width, height]. |
| `export.plot_unpaired` | `false` | When `true`, single-spectrum plots are generated for queries with no reference. When `false`, those queries are skipped. |

### `performance:`

| Key | Default | Description |
|---|---|---|
| `performance.workers` | `4` | Parallel workers. Same knob as `settings.json > performance.workers`. |

### `koina:`

| Key | Default | Description |
|---|---|---|
| `koina.enabled` | `true` | Master switch. `--no-koina` flips it off. |
| `koina.model` | `"ms2pip_timsTOF2024"` | Model ID. |
| `koina.collision_energy` | `35` | NCE. |
| `koina.instrument_type` | `"TIMSTOF"` | Instrument hint. |

---

## Command-line flags

Every flag overrides the corresponding YAML key for a single invocation.

| Flag | Argument | Overrides |
|---|---|---|
| `-c`, `--config` | PATH | *(required)* path to YAML config |
| `--idxml` | PATTERN… | `input.idxml_files` |
| `--mzml` | PATTERN… | `input.mzml_files` |
| `--peptides` | PATH_OR_ALL | Top-level `peptides:` |
| `--reference-substring` | STR | `reference.substring` |
| `--ion-types` | ION… | `annotation.ion_types` |
| `--mz-tolerance` | DA | `annotation.mz_tolerance` |
| `--export-dir` | DIR | `export.directory` |
| `--project-name` | NAME | `export.project_name` |
| `--formats` | FMT… | `export.formats` (choices: png svg pdf html) |
| `--workers` | N | `performance.workers` |
| `--no-plots` | *(flag)* | Skip plot generation; export CSV/Parquet only |
| `--no-koina` | *(flag)* | Disable Koina regardless of YAML |
| `--koina-model` | ID | `koina.model` |
| `--koina-ce` | EV | `koina.collision_energy` |
| `--koina-instrument` | TYPE | `koina.instrument_type` |
| `-v`, `--verbose` | *(flag)* | Enable DEBUG logging |
| `-q`, `--quiet` | *(flag)* | Suppress non-error output |

---

## See also

- [../cli/usage.md](../cli/usage.md) — invocation, exit codes, recipes.
- [app.md](app.md) — the `settings.json` layer that seeds CLI defaults.
- [overview.md](overview.md) — three-layer precedence rules.
