# CLI Recipes & Debugging

Common invocation patterns and how to recover when something fails.

For the full flag reference, see [../settings/cli.md](../settings/cli.md).

---

## Peptide specification

The `peptides:` YAML key (or `--peptides` flag) accepts three forms:

| Form | Example | Behaviour |
|---|---|---|
| Keyword | `peptides: all` | Scan every idXML and use every unique peptide found, capped at 10 000. Progress bar appears during the scan. |
| File path | `peptides: /data/list.txt` | Read from disk. Comma- or newline-separated. Case is folded to uppercase. |
| Inline list | `peptides: [ASDDFDQASY, RTEFNLNQY]` | YAML flow-style list. |

All three route through `cli.config.load_peptides()`. The resolved list
and mode are recorded in `run_settings.json` for reproducibility.

---

## Reference substring and single-spectrum mode

The `reference.substring` value splits your selected mzML files into two
groups:

- **Reference files** — name contains the substring. Used to build the
  reference-intensity side of each mirror plot.
- **Query files** — everything else. Used as the top side of each mirror
  plot.

Leave the substring empty to enter **single-spectrum mode**: queries are
plotted alone (no mirror), no pairwise metrics are computed, and no
Koina fallback is attempted. Useful for one-off inspection of a batch.

When Koina is enabled *and* a query has no matching experimental
reference, MSpecViz automatically predicts a reference from the peptide's
amino-acid sequence. The predicted spectrum is labelled `Koina` in the
output.

---

## Common recipes

### Re-run with a tweaked tolerance

```bash
python -m cli -c my_config.yaml --mz-tolerance 0.05
```

Creates a new numbered export directory because the settings fingerprint
changed.

### Sweep peptide lists from a shell loop

```bash
for PEPS in lists/*.txt; do
  python -m cli -c base_config.yaml \
    --peptides "$PEPS" \
    --project-name "sweep_$(basename "$PEPS" .txt)" \
    -q
done
```

### Regenerate HTML plots from an already-processed dataset

Make a copy of the config, change `export.formats` to `[html]`, and re-run:

```yaml
export:
  formats: [html]
```

```bash
python -m cli -c html_only.yaml
```

The DataFrames are rebuilt but plot generation is limited to HTML.

### Pipeline that touches disk but not Koina

```bash
python -m cli -c my_config.yaml --no-koina
```

Any query with no experimental reference is skipped (or emitted as a
single-spectrum plot when `export.plot_unpaired: true`).

### Headless validation on CI

```bash
python -m cli -c cli/config/config_example.yaml --no-plots -q
```

Exits with 0 on success; 1 on pipeline error; 2 on config error. Safe
for CI smoke tests — no GUI, no plot rendering.

---

## Debugging tips

- **`--workers 1`** turns off the `ProcessPoolExecutor` and runs the
  pipeline sequentially. Crashes surface as plain tracebacks instead of
  silent worker exits.
- **`-v`** prints DEBUG-level logs from `pyopenms`, `koinapy`, and every
  `mspecviz.*` module. Useful when a Koina call silently returns no
  spectrum.
- **Stale cache** — if you've fixed a loading bug and re-runs keep using
  the buggy output, delete the workspace's `.mspecviz_cache/` directory.
  The CLI doesn't invalidate cache automatically. See
  [../settings/app.md#cache](../settings/app.md#cache).
- **Same-named files** in different directories are disambiguated with a
  parent-directory suffix by `mspecviz.path_utils`. Look for names like
  `foo.mzML [experimentA]` in the output if you notice them.
- **`run_settings.json`** is the canonical record of what actually ran.
  If a result looks wrong, start there — the resolved peptide list,
  expanded glob paths, and the final merged settings are all captured
  verbatim.
