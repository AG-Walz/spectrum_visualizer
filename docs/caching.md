# Parquet Caching System

MSpecViz automatically caches loaded spectra in Parquet format so that repeated runs with the same inputs load ~100x faster (skipping all mzML XML parsing and index building).

## How It Works

### Fingerprinting

Before loading spectra, MSpecViz computes a **fingerprint** — a 12-character hex string (e.g., `48020b32ce9d`) derived from a SHA-256 hash of:

1. **All file paths** (idXML + mzML), sorted, each combined with the file's modification time (`mtime`) and size on disk
2. **All target peptide sequences**, sorted and uppercased
3. **The reference substring**, lowercased

If **any** of these change — different files selected, a peptide added/removed, an mzML file reprocessed (changing its mtime), or a different reference substring — the fingerprint changes and the cache is not reused.

> **Note:** The fingerprint uses filesystem metadata (mtime + size), not file content. This is fast (no need to read multi-GB files) but means that copying a file to a new path or touching it without changing content will produce a different fingerprint.

### Cache Lookup

The app checks for a file named `{fingerprint}.parquet` in the cache directory. If found and non-empty, it loads from Parquet instead of mzML.

### Cache Storage

Cache files are stored in:

```
{export_dir}/.mspecviz_cache/
├── 48020b32ce9d.parquet    # Cached spectra (binary, columnar)
└── 48020b32ce9d.json       # Human-readable sidecar metadata
```

The default location is `exports/.mspecviz_cache/`. This is configurable via `settings.json`:

```json
{
  "cache": {
    "enabled": true,
    "subdirectory": ".mspecviz_cache"
  }
}
```

### What's Stored

Each row in the Parquet file represents one loaded spectrum:

| Column | Type | Description |
|--------|------|-------------|
| `source_mzml` | string | Source mzML display name (disambiguated) |
| `spectrum_id` | string | Native spectrum ID (e.g., `index=12345`) |
| `peptidoform` | string | Full sequence with modifications |
| `base_sequence` | string | Unmodified peptide sequence |
| `variant_key` | string | Grouping key (`sequence\|modifications`) |
| `charge` | int | Precursor charge state |
| `precursor_mz` | float | Precursor m/z |
| `rt` | float | Retention time |
| `is_reference` | bool | Whether this is a reference spectrum |
| `mz` | float[] | m/z array |
| `intensity` | float[] | Intensity array |
| `score` | float | PSM score |
| `idxml_file` | string | Source idXML display name |

The sidecar `.json` file contains metadata for inspection:

```json
{
  "created_by": "MSpecViz auto-cache",
  "created_at": "2026-02-13T13:46:20",
  "total_spectra": 11,
  "reference_spectra": 5,
  "query_spectra": 6,
  "fingerprint": "48020b32ce9d",
  "idxml_files": ["..."],
  "mzml_files": ["..."],
  "target_peptides": ["PEPTIDEK", "..."],
  "reference_substring": "synp"
}
```

## When the Cache Is Saved

The cache is written **only after a complete, successful run**. The save happens at the very end of the "Process & Export" step, after all spectra have been loaded, annotated, and exported.

## Interrupted Runs

**If you stop a run mid-process, no cache is saved.** The spectra loaded up to that point exist only in memory and are discarded. The next run will start over from scratch, loading all spectra from mzML again.

This is by design — a partial cache could contain incomplete data (e.g., missing spectra from files that weren't processed yet), which would silently produce incorrect results on subsequent runs.

## Cache Invalidation

The cache is automatically invalidated (not reused) when:

- A different set of idXML or mzML files is selected
- The peptide list changes
- The reference substring changes
- Any input file is modified on disk (changing its mtime or size)

To **manually invalidate** the cache, delete the `.parquet` and `.json` files from the cache directory:

```bash
rm exports/.mspecviz_cache/*.parquet exports/.mspecviz_cache/*.json
```

## Disabling the Cache

Set `cache.enabled` to `false` in `settings.json`:

```json
{
  "cache": {
    "enabled": false
  }
}
```

## Troubleshooting

**Cache shows stale results after a code change:**
If you fix a bug in spectrum loading or annotation, the cache will still contain data from the buggy code (the fingerprint only tracks input files, not code). Delete the cache files after any code change that affects spectrum loading:

```bash
rm exports/.mspecviz_cache/*.parquet exports/.mspecviz_cache/*.json
```

**Cache not being created:**
- Check that the export directory is writable
- Check the Streamlit logs for "Cache save failed" warnings
- Ensure the run completed fully (interrupted runs don't save)
