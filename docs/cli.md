# CLI

`python -m cli` runs the full pipeline headlessly from a YAML config, writing the
same versioned export directory as the Streamlit Export tab (see
[Export → directory layout and versioning](export.md#directory-layout-and-versioning)).

---

## Invocation

```bash
python -m cli cli/config/config_example.yaml     # from the project root
```

One positional argument: the YAML config path. There are **no** command-line flag
overrides — copy the YAML and edit it for per-run tweaks. Log level is set out of
band:

```bash
MSPECVIZ_LOG_LEVEL=ERROR python -m cli my_config.yaml   # quiet — errors only
MSPECVIZ_LOG_LEVEL=DEBUG python -m cli my_config.yaml   # full DEBUG from pyopenms/koinapy/core
```

Each stage logs a header and a `tqdm` bar (driven by the same `run_tasks`
callbacks the Streamlit progress bars use). Exit code is `0` on success and `1` on
any failure (`sys.exit("Error: …")` — there is no distinct config-error code; the
test named `test_missing_config_exits_two` only asserts non-zero).

---

## YAML config

Flat `snake_case` keys, each matching a `DEFAULTS` field exactly (precedence:
`DEFAULTS < YAML`; see [Settings → sources and precedence](settings.md#sources-and-precedence)).
Required: `idxml_files`, `mzml_files`, `peptide_source`;
everything else inherits from `DEFAULTS`. The bundled
`cli/config/config_example.yaml` is a fully-annotated template listing every
field — **copy it rather than hand-writing one**. Minimum viable config:

```yaml
idxml_files:
  - /path/to/results.idXML
mzml_files:
  - /path/to/raw/*.mzML          # globs expanded with recursive=True (** works)
peptide_source: /path/to/peptides.txt   # file (one peptide/line) or inline [PEP1, PEP2]
```

Notes that don't fit in inline comments:

- `summary_metric` must be a key in `core/plotting/summary.py:METRIC_COLUMNS`.
- `summary_plot_types` accepts a subset of `{dual_axis, ridgeline}` — no other type
  is implemented (see [Export → plots](export.md#plots)).
- `peptide_source` peptides are uppercased and whitespace-stripped before lookup.

### Reference substring and single-spectrum mode

`reference_substring` splits the selected mzMLs by a **case-sensitive** `in` test on
the basename: matches are reference files (bottom of each mirror plot), the rest are
queries. Pick something distinctive (`Synpep`, `_Ref_`) — a short substring like `1`
matches almost anything. Leave it empty for **single-spectrum mode**: queries are
plotted alone, no pairwise metrics, no Koina fallback. When Koina is enabled and a
query has no experimental reference, a reference is predicted and labelled `Koina`.

---

## Recipes

- **Tweak a parameter** (e.g. `mz_tolerance: 0.05`): edit a copy and re-run. A
  changed setting fingerprint allocates a new numbered export dir (see
  [Export → directory layout and versioning](export.md#directory-layout-and-versioning)).
- **HTML only:** set `export_formats: [html]` — DataFrames rebuild, plots limited to HTML.
- **No Koina:** `koina_enabled: false` — unpaired queries are skipped unless
  `plot_unpaired_queries: true`.
- **CI smoke test:** `no_plots: true` with `MSPECVIZ_LOG_LEVEL=ERROR` — exports data
  only, exits non-zero on any error.
- **Sweep peptide lists:** generate one YAML per list in a shell loop (`sed` the
  `peptide_source` / `project_name` lines), run each.

---

## Debugging

- **`workers: 1`** forces the sequential path in `run_tasks` — crashes surface as
  real tracebacks instead of a generic "worker crashed".
- **`MSPECVIZ_LOG_LEVEL=DEBUG`** when a Koina call silently returns no spectrum.
- **Stale CLI plot** is never a cache issue (the CLI doesn't cache) — it's the
  versioned-export fingerprint reusing a directory; change a parameter to force a
  new one.
- **Same-named files** across directories: the CLI does *not* call
  `ensure_unique_basenames`, so colliding globs silently overwrite — review your
  globs when sweeping directories.
- **`run_settings.json`** is the canonical record of what ran (resolved peptides,
  expanded globs, merged settings) — start there when a result looks wrong.
