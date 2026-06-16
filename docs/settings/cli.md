# CLI Settings — YAML

The single configuration surface for `python -m cli` is a YAML file
(passed as the only positional argument). There are no command-line
flag overrides — copy the YAML and edit per-run instead.

Every YAML key matches a field on the merged runtime config exactly —
no translation table. Defaults come from
`core.defaults.DEFAULTS`. See [overview.md](overview.md) for
precedence; [../cli/usage.md](../cli/usage.md) for invocation.

The log level is set via the `MSPECVIZ_LOG_LEVEL` environment variable
(default `INFO`).

---

## YAML config structure

Flat snake_case keys. Every key is optional except `idxml_files`,
`mzml_files`, and `peptide_source`.

```yaml
# Inputs
idxml_files: []                       # list of paths or glob patterns
mzml_files:  []                       # list of paths or glob patterns
peptide_source: /path/to/peptides.txt  # file path | [PEP1, PEP2]
reference_substring: ""               # empty = single-spectrum mode

# Annotation
ion_types: [b, y]                     # any of: a b c x y z precursor immonium
mz_tolerance: 0.02                    # Da
add_neutral_losses: false             # H2O + NH3 losses
intensity_mode: Relative              # "Relative" or "Raw"

# Export
export_directory: ./exports           # base output dir
project_name: MSpecViz_analysis       # subdir prefix
export_formats: [png]                 # png svg pdf html (any combo)
export_resolution: 300                # raster DPI
figure_size: [12, 8]                  # inches [width, height]
plot_unpaired_queries: false          # emit single-spectrum plots when no ref
no_plots: false                       # skip plot generation; export only data

# Performance
workers: 4                            # 1 = sequential fallback

# Summary plots
summary_metric: "R² binned"           # default summary-plot metric
summary_plot_types: [dual_axis, ridgeline]  # any subset of {dual_axis, ridgeline}

# Koina AI prediction
koina_enabled: false                  # master switch for AI prediction
koina_model: ms2pip_timsTOF2024       # any id in core/defaults.py:KOINA_MODELS
koina_collision_energy: 35            # NCE
koina_instrument_type: TIMSTOF        # QE | LUMOS | TIMSTOF | SCIEXTOF
```

A few notes that don't fit in inline comments:

- `idxml_files` / `mzml_files` globs are expanded with `recursive=True`.
- `summary_metric` must be a key in `core/plotting/summary.py:METRIC_COLUMNS`.
- `summary_plot_types` accepts a subset of `{dual_axis, ridgeline}`; no
  other type is implemented.

---

### `peptide_source` interpretation

| Value | Behaviour |
|---|---|
| `/abs/path/peptides.txt` | Peptide list read from disk — one peptide per line. |
| `[PEPTIDE1, PEPTIDE2]` | Inline list (YAML flow style). |

Peptides are case-folded to uppercase and stripped of whitespace before
the lookup.

---

## Environment variables

| Variable | Default | Description |
|---|---|---|
| `MSPECVIZ_LOG_LEVEL` | `INFO` | Python logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. Lower-case values are uppercased automatically. |

---

## See also

- [../cli/usage.md](../cli/usage.md) — invocation, exit codes, recipes.
- [app.md](app.md) — the `settings.json` layer (Streamlit-only).
- [overview.md](overview.md) — precedence rules.
