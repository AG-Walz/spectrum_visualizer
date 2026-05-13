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

- **Reference files** — basename contains the substring (case-sensitive
  Python `in` test on the filename, not the full path).
- **Query files** — everything else. Used as the top side of each mirror
  plot.

Pick a substring distinctive enough to be unambiguous (e.g. `Synpep`,
`Standard`, `_Ref_`). A short substring like `1` will match almost
anything; trim or pre-rename your files if needed.

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

Exits with 0 on success and non-zero (always 1) on any error. Safe for
CI smoke tests — no GUI, no plot rendering.

---

## Debugging tips

- **`workers: 1`** in YAML turns off the `ProcessPoolExecutor` and runs
  the pipeline sequentially. Crashes surface as plain tracebacks instead
  of silent worker exits.
- **`MSPECVIZ_LOG_LEVEL=DEBUG`** prints DEBUG-level logs from `pyopenms`,
  `koinapy`, and every `core.*` module. Useful when a Koina call
  silently returns no spectrum.
- **Stale cache** — only relevant to the Streamlit app, which persists
  loaded spectra to `<workspace>/.mspecviz_cache/`. The CLI does **not**
  cache and always re-runs the load stage. If you see a stale plot from
  the CLI, the cause is the versioned-export fingerprint reusing an
  existing directory — delete that directory or change a parameter to
  force a new fingerprint.
- **Same-named files** — if you select two files with the same basename
  from different directories, both the Streamlit Setup tab and any code
  path that calls `core.common.helpers.ensure_unique_basenames`
  will raise `DuplicateBasenameError`. Rename or stage the files
  separately. The CLI does not currently call this check itself, so
  glob results with collisions silently overwrite — review your globs
  if you're sweeping across directories.
- **`run_settings.json`** is the canonical record of what actually ran.
  If a result looks wrong, start there — the resolved peptide list,
  expanded glob paths, and the final merged settings are all captured
  verbatim. Reproducing a Streamlit run from the CLI is just a matter of
  translating that JSON into YAML.
