# File Selection

This document explains how MSpecViz handles file selection, directory browsing, and disambiguation of files with identical names.

---

## Supported File Types

### idXML — Identification Files

idXML is the OpenMS XML format for peptide identifications. Each file contains:

- **PeptideIdentification** records with one or more **PeptideHit** entries
- Each PeptideHit has: amino acid sequence, modifications, charge state, score
- **SpectraData** metadata mapping spectrum references to source mzML files via `id_merge_index`

idXML files are produced by OpenMS search engines (MSGF+, Comet, etc.) or converted from other formats (mzIdentML, pepXML) using OpenMS tools like `IDFileConverter`.

### mzML — Spectrum Files

mzML is the PSI standard format for mass spectrometry data. Each file contains:

- Thousands to millions of MS/MS spectra
- Each spectrum stores m/z and intensity arrays in compressed binary format
- An index at the end of the file enables fast random access to individual spectra
- **Native IDs** identify each spectrum (e.g., `index=12345`, `scan=6789`)

mzML files are typically 1–10 GB. MSpecViz accesses them via indexed random access — only the spectra you need are loaded, not the entire file.

---

## Directory Browser

Each file type has its own directory browser. The browser lets you navigate the server filesystem and select files without uploading.

### Navigation

- **Path input**: Type or paste a directory path. Press Enter to navigate.
- **Up button (⬆️)**: Go to the parent directory.
- **Folder grid**: Click a folder name to navigate into it. Up to 15 folders are shown; if more exist, a count indicates how many are hidden.

### File Selection

1. The browser lists all files matching the expected extension (`.idXML` or `.mzML`)
2. Use the **filter box** to narrow the list by typing part of a filename
3. Use the **Select All** checkbox to select every file in the current directory
4. Or use the **multi-select dropdown** to pick individual files

### Accumulating Across Directories

Selections persist as you navigate between directories. If you select 3 files in directory A, then navigate to directory B and select 2 more, you'll have 5 files total.

The **"All selected files"** expander (at the bottom of each browser) shows all selections grouped by directory. You can:

- **Remove individual files** with the ❌ button
- **Clear all selections** with the 🗑️ button

---

## Path Disambiguation

When you select files with identical filenames from different directories, MSpecViz needs to tell them apart. It uses a **minimum-depth disambiguation algorithm**:

### Algorithm

1. Start with just the filename for each path
2. Check if any names collide (same string for different paths)
3. For colliding names, add one parent directory component
4. Repeat until all names are unique

### Example

Given these selected files:

```
/data/batch1/run1/sample.mzML
/data/batch1/run2/sample.mzML
/data/batch2/run1/sample.mzML
/data/other/unique_file.mzML
```

The disambiguation produces:

```
sample.mzML           → still collides (3 files)
run1/sample.mzML      → still collides (2 files share this)
batch1/run1/sample.mzML → unique ✓
batch1/run2/sample.mzML → unique ✓ (already resolved at depth 2)
batch2/run1/sample.mzML → unique ✓ (already resolved at depth 2)
unique_file.mzML       → unique from the start ✓
```

Final display names:

| Full Path | Display Name |
|-----------|-------------|
| `/data/batch1/run1/sample.mzML` | `batch1/run1/sample.mzML` |
| `/data/batch1/run2/sample.mzML` | `run2/sample.mzML` |
| `/data/batch2/run1/sample.mzML` | `batch2/run1/sample.mzML` |
| `/data/other/unique_file.mzML` | `unique_file.mzML` |

### Deduplication

If the same file is selected twice (exact same path), only one copy is kept. This happens silently — no warning is shown.

---

## Parquet Shortcut

Below both file browsers, there is a **Pre-converted Parquet** section. If you have a `.parquet` file from a previous analysis run (auto-cached or manually converted), paste its full path here.

When a valid Parquet file is selected:

- Steps 2–5 (Peptides, Reference, Settings, Search) become optional
- You can go directly to Step 6 (Process)
- The Parquet file already contains all spectra, peptide information, and reference/query categorization

This is the fastest way to rerun an analysis — loading takes under 1 second regardless of dataset size.

### How to Get a Parquet File

1. **Automatic**: Complete a full run once. MSpecViz saves a cache in `exports/.mspecviz_cache/`. Copy the `.parquet` file to use it.
2. **In-app**: After Step 5 (Search), click "Convert to Parquet".
3. **CLI tool**: Use `python tools/convert_mzml_to_parquet.py` for batch conversion. See [tools/README.md](../tools/README.md).

---

## mzML Resolution

A key challenge is linking each PSM in an idXML file to the correct mzML file. The idXML stores spectrum references like `index=12345`, but this index is meaningless without knowing *which* mzML file it refers to.

### How Resolution Works

1. The idXML's `ProteinIdentification` section contains `SpectraData` entries listing the source mzML files
2. Each `PeptideIdentification` has an `id_merge_index` pointing to one of these entries
3. MSpecViz extracts the mzML filename from the `SpectraData` entry and matches it against your selected mzML files

### When Multiple Files Share a Name

If the metadata says `sample.mzML` but you have three files named `sample.mzML` from different directories, MSpecViz uses these heuristics (in order):

1. **Path suffix matching**: Compare the full metadata path against each candidate. The best suffix match wins.
2. **idXML proximity**: Prefer the mzML file in the same directory (or closest directory) to the idXML that produced the PSM.
3. **First candidate fallback**: If all else fails, use the first matching file (with a warning).

---

## Validation

Step 1 requires one of:

| Condition | Outcome |
|-----------|---------|
| ≥1 idXML + ≥1 mzML selected | Proceed to Step 2 |
| Valid Parquet file path | Skip to Step 6 |
| Neither | Warning message, cannot proceed |
