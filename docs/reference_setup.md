# Reference Setup

This document explains how MSpecViz categorizes mzML files into **reference** and **query** groups, and the implications for visualization.

---

## Concept

In a typical MSpecViz workflow, you compare experimental spectra against a known standard:

- **Reference spectra**: From a trusted source — synthetic peptide library, well-characterized sample, or a specific experimental condition you want to compare against
- **Query spectra**: The spectra you want to evaluate — experimental samples, test conditions, replicates

Mirror plots place the query spectrum on **top** and the reference on the **bottom**, making it easy to visually assess agreement.

---

## How Categorization Works

MSpecViz uses a **substring match** on mzML filenames to categorize files:

1. You enter a reference identifier string (e.g., `Synpep`)
2. For each selected mzML file, MSpecViz checks if the filename contains this string (case-insensitive)
3. Matching files → **Reference**
4. Non-matching files → **Query**

### Example

Selected files:
```
Exp01_Synpep_Rep01.mzML     → contains "Synpep" → Reference
Exp01_Synpep_Rep02.mzML     → contains "Synpep" → Reference
Exp01_Sample_A.mzML         → does not contain → Query
Exp01_Sample_B.mzML         → does not contain → Query
Exp01_Sample_C.mzML         → does not contain → Query
```

### Case Insensitivity

The match is case-insensitive: entering `synpep`, `Synpep`, or `SYNPEP` all produce the same result.

---

## Preview

After entering a substring, a two-column preview displays:

- **Left column**: Reference files (with count)
- **Right column**: Query files (with count)
- Up to 5 files are shown per column; if more exist, a count indicates the remainder

This lets you verify the categorization before proceeding.

---

## What Happens During Processing

When spectra are loaded (Step 6), the reference substring determines processing order and comparison pairing:

1. **Reference spectra are loaded first**: For each peptide variant, all reference spectra from all reference mzML files are loaded and annotated. These are cached in memory.

2. **Query spectra are loaded one at a time**: Each query spectrum is compared against **every** reference spectrum for the same peptide variant. This means:
   - 1 query + 3 references = 3 mirror plots
   - 5 queries + 3 references = 15 mirror plots
   - 10 queries + 2 references = 20 mirror plots

3. **Matching is by peptide variant**: A query spectrum for `PEPTIDEK|M3[Ox]` is only compared to reference spectra with the same variant key (`PEPTIDEK|M3[Ox]`), not to `PEPTIDEK` without oxidation.

---

## No Reference Substring

If you **leave the reference field blank**, MSpecViz operates in **individual annotation mode**:

- No mirror plots are generated
- Each spectrum is plotted individually with theoretical fragment annotations
- All spectra are treated as queries
- Useful for simply visualizing fragmentation patterns without comparison

---

## Unpaired Spectra

A "pairing" happens when a peptide variant has at least one reference **and** at least one query spectrum. Common situations where pairing fails:

| Situation | What Happens |
|-----------|-------------|
| Peptide found in query files but **not in any reference file** | Query spectra have no reference to compare against |
| Peptide found in reference files but **not in any query file** | Reference spectra are loaded but never used |
| Peptide found in both | Normal mirror plot generation |

### Plot Unpaired Option

The **"Plot unpaired spectra"** checkbox (in Step 6 Export Settings) controls what happens when a query has no reference:

- **Off (default)**: The query spectrum is silently skipped. Only mirror plots (with both query and reference) are exported.
- **On**: A single annotated spectrum plot is created for the unpaired query. No comparison metrics are computed.

---

## Choosing a Good Reference Substring

Tips for choosing an effective substring:

1. **Be specific enough** to avoid false matches. If query files also contain "ref" in their names, use a more distinctive string.

2. **Be general enough** to match all reference files. If references are named `SynPep_01.mzML` through `SynPep_05.mzML`, entering `Synp` is sufficient.

3. **Check the preview** before proceeding. The two-column display immediately shows whether your substring produces the expected split.

4. **Common patterns**:
   - Synthetic peptide libraries: `Synpep`, `SynPep`, `synthetic`
   - Standard references: `Std`, `Standard`, `Ref`
   - Control samples: `Control`, `Ctrl`, `Blank`
   - Specific conditions: `Treatment`, `Drug`, `HighDose`

---

## Parquet Files and References

When loading from a Parquet file, the reference categorization is already embedded in the data (the `is_reference` column). However, if you specify a reference substring in Step 3, MSpecViz will use it to re-categorize based on the `source_mzml` column in the Parquet file. This lets you reanalyze the same data with a different reference definition without re-extracting spectra.
