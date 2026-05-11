# CLI Recipes & Debugging

Common invocation patterns and how to recover when something fails.

For the full YAML schema, see [../settings/cli.md](../settings/cli.md).

---

## Peptide specification

The `peptide_source:` YAML key accepts two forms:

| Form | Example | Behaviour |
| --- | --- | --- |
| File path | `peptide_source: /data/list.txt` | Read from disk. Comma- or newline-separated. Case is folded to uppercase. |
| Inline list | `peptide_source: [ASDDFDQASY, RTEFNLNQY]` | YAML flow-style list. |

Both route through `cli.main._load_peptides()`. The resolved list is
recorded in `run_settings.json` for reproducibility.

---

## Reference substring and single-spectrum mode

The `reference_substring` value splits your selected mzML files into two
groups:

- **Reference files** — name contains the substring. Used to build the
  reference-intensity side of each mirror plot.
- **Query files** — everything else. Used as the top side of each mirror
  plot.

Leave the substring empty to enter **single-spectrum mode**: queries are
plotted alone (no mirror), no pairwise metrics are computed, and no
Koina fallback is attempted. Useful for one-off inspection of a batch.

When Koina is enabled *and* a query has no matching experimental
reference, MSpecViz automatically predicts a reference from the
peptide's amino-acid sequence. The predicted spectrum is labelled
`Koina` in the output.

---

## Common recipes

### Re-run with a tweaked tolerance

Copy your YAML and edit `mz_tolerance`:

```yaml
mz_tolerance: 0.05
```

```bash
python -m cli my_tweaked_config.yaml
```

Creates a new numbered export directory because the settings
fingerprint changed.

### Sweep peptide lists from a shell loop

Generate a temporary YAML per list:

```bash
for PEPS in lists/*.txt; do
  NAME=$(basename "$PEPS" .txt)
  sed -e "s|peptide_source:.*|peptide_source: $PEPS|" \
      -e "s|project_name:.*|project_name: sweep_$NAME|" \
      base_config.yaml > "/tmp/cfg_$NAME.yaml"
  MSPECVIZ_LOG_LEVEL=ERROR python -m cli "/tmp/cfg_$NAME.yaml"
done
```

### Regenerate HTML plots from an already-processed dataset

Make a copy of the config, change `export_formats` to `[html]`, and re-run:

```yaml
export_formats: [html]
```

```bash
python -m cli html_only.yaml
```

The DataFrames are rebuilt but plot generation is limited to HTML.

### Pipeline that touches disk but not Koina

Set `koina_enabled: false` in YAML:

```yaml
koina_enabled: false
```

Any query with no experimental reference is skipped (or emitted as a
single-spectrum plot when `plot_unpaired_queries: true`).

### Headless validation on CI

```yaml
no_plots: true
```

```bash
MSPECVIZ_LOG_LEVEL=ERROR python -m cli cli/config/config_example.yaml
```

Exits with 0 on success, 1 on pipeline error, 2 on config error. Safe
for CI smoke tests — no GUI, no plot rendering.

---

## Debugging tips

- **`workers: 1`** in YAML turns off the `ProcessPoolExecutor` and runs
  the pipeline sequentially. Crashes surface as plain tracebacks instead
  of silent worker exits.
- **`MSPECVIZ_LOG_LEVEL=DEBUG`** prints DEBUG-level logs from `pyopenms`,
  `koinapy`, and every `mspecviz.*` module. Useful when a Koina call
  silently returns no spectrum.
- **Stale cache** — if you've fixed a loading bug and re-runs keep using
  the buggy output, delete the workspace's `.mspecviz_cache/` directory.
  The CLI doesn't invalidate cache automatically. See
  [../settings/app.md#cache](../settings/app.md#cache).
- **Same-named files** in different directories are disambiguated with a
  parent-directory suffix by `mspecviz.common.path_utils`. Look for
  names like `foo.mzML [experimentA]` in the output if you notice them.
- **`run_settings.json`** is the canonical record of what actually ran.
  If a result looks wrong, start there — the resolved peptide list,
  expanded glob paths, and the final merged settings are all captured
  verbatim.
