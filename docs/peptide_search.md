# Peptide Search

This document explains how MSpecViz searches idXML files for target peptides and resolves the corresponding mzML spectra.

---

## Overview

The peptide search is the bridge between your input files and the analysis. It:

1. Parses each idXML file for peptide-spectrum matches (PSMs)
2. Groups PSMs by peptide variant (base sequence + modification combination)
3. Resolves which mzML file contains each spectrum
4. Returns a mapping of variant keys to match lists

---

## Input: Target Peptides

You provide a list of **base sequences** (unmodified amino acid sequences). MSpecViz searches for these sequences and automatically finds all modified variants.

For example, if you enter `PEPTIDEK`, the search will find:
- `PEPTIDEK` (unmodified)
- `PEPTIDEK` with Oxidation at position 3
- `PEPTIDEK` with Phosphorylation at position 5
- Any other modification variant present in the idXML

**Case insensitive**: `peptidek`, `PEPTIDEK`, and `Peptidek` are all equivalent.

---

## How the Search Works

### Step 1: Load idXML

Each idXML file is loaded using PyOpenMS. The loader extracts:
- All `PeptideIdentification` records
- The `ProteinIdentification` headers (which contain `SpectraData` metadata)

### Step 2: Extract SpectraData Map

The `SpectraData` entries in the idXML map integer indices to mzML file information:

```
id_merge_index 0  →  /path/to/sample1.mzML
id_merge_index 1  →  /path/to/sample2.mzML
id_merge_index 2  →  /path/to/Synpep_ref.mzML
```

This mapping is essential for resolving which mzML file contains each spectrum.

### Step 3: Match Peptides

For each `PeptideIdentification`, the search:

1. Extracts the spectrum native ID (e.g., `index=12345`)
2. Looks up the `id_merge_index` to determine the source mzML
3. Gets the top-ranked `PeptideHit`
4. Extracts the base sequence and checks if it matches any target peptide
5. If matched, extracts: peptidoform (full modified sequence), charge, score, modifications

### Step 4: Group by Variant

Matches are grouped into **VariantGroups** — one group per unique combination of base sequence and modification pattern:

```
PEPTIDEK|                → unmodified variant
PEPTIDEK|M3[Ox]          → oxidized variant
```

Each VariantGroup contains a list of PeptideMatch objects with all the details needed to load the corresponding spectrum.

---

## Modification Extraction

Modifications are extracted from each PeptideHit and formatted as compact inline strings using abbreviations for common modifications:

| Modification Location | Format Example |
|-----------------------|---------------|
| N-terminal | `NTerm[Ac]` |
| Residue-specific | `M5[Ox]` (Oxidation at position 5) |
| C-terminal | `CTerm[Amidated]` |

Multiple modifications are joined with `|`:
```
M5[Ox]|S8[Ph]
```

Common abbreviations: Ox (Oxidation), Ph (Phospho), CAM (Carbamidomethyl), Ac (Acetyl), Deam (Deamidated), Me (Methylation). Unknown modifications use their full name.

The **display label** shows modifications inline within the sequence: `PEPTM[Ox]DEK` instead of `PEPTMDEK [Oxidation@M5]`.

The **modification key** (used for grouping) normalizes modifications so that the same chemical modification at the same position always produces the same string, regardless of idXML formatting differences.

---

## mzML Resolution

After finding PSMs, each must be mapped to a specific mzML file on disk. This is non-trivial when files share names.

### Simple Case: Unique Filenames

If the metadata says `sample1.mzML` and only one of your selected files has that name, the mapping is trivial.

### Complex Case: Duplicate Filenames

When multiple selected mzML files share the same filename, MSpecViz uses a three-strategy resolution:

1. **Path suffix matching**: The idXML metadata often contains a longer path (e.g., `/data/experiment1/sample.mzML`). MSpecViz compares this against each candidate's full path and picks the best suffix match.

2. **idXML directory proximity**: If suffix matching is inconclusive, prefer the mzML file closest (in directory structure) to the idXML that produced the PSM.

3. **First candidate fallback**: If all strategies are inconclusive, use the first matching file and log a warning.

### PeptideMatch Fields

Each resolved match contains:

| Field | Description |
|-------|-------------|
| `peptidoform` | Full modified sequence (e.g., `PEPTM(Oxidation)IDEK`) |
| `base_sequence` | Unmodified sequence (e.g., `PEPTMIDEK`) |
| `spectrum_id` | Native spectrum ID (e.g., `index=12345`) |
| `charge` | Precursor charge state |
| `score` | PSM confidence score |
| `mzml_filename` | Bare mzML filename from metadata |
| `mzml_resolved_path` | Full resolved path to the mzML file on disk |
| `idxml_file` | Disambiguated idXML display name |
| `idxml_resolved_path` | Full path to the source idXML file |

---

## Results Table

The search results are summarized in a table with one row per peptide variant:

| Column | Meaning |
|--------|---------|
| **Peptide** | Base amino acid sequence |
| **Modification** | Modification string (empty if unmodified) |
| **PSMs** | Total number of PSMs found for this variant |
| **idXML Files** | How many idXML files contained this variant |
| **mzML Resolved** | How many PSMs were successfully mapped to mzML files |
| **Unresolved** | PSMs where the mzML could not be determined |

### Unresolved PSMs

A PSM is "unresolved" when:
- The `id_merge_index` is missing from the `PeptideIdentification`
- The `SpectraData` list in the idXML is empty
- The mzML filename from metadata doesn't match any selected file

Unresolved PSMs are skipped during processing. If many PSMs are unresolved, verify that all relevant mzML files are selected.

---

## Convert to Parquet

After a successful search, you can convert the matched spectra to Parquet format for faster loading:

1. A text input shows a default output path (in the export directory)
2. Click **Convert to Parquet** to start conversion
3. Progress is shown as a caption below the button
4. When complete, the Parquet file path is stored in session state

The conversion:
- Uses `OnDiscMSExperiment` for memory-efficient mzML access
- Extracts only the matched spectra (not the entire file)
- Writes a compressed Parquet file with all metadata
- Can be used in Step 1 on future runs for instant loading

See [Parquet Conversion Tool](../tools/README.md) for the standalone CLI version.
