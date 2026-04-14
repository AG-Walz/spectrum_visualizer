# Metrics Reference

This document describes all quantitative metrics computed by MSpecViz when
comparing query and reference spectra. Implementations are in
`mspecviz/spectrum_metrics.py` (7 comparison metrics) and
`mspecviz/plot_worker.py` (17 ion match statistics).

---

## Table of Contents

1. [Overview](#overview)
2. [7 Comparison Metrics](#7-comparison-metrics)
   - [R-squared Matching (R2m)](#r-squared-matching-r2m)
   - [R-squared Binned (R2b)](#r-squared-binned-r2b)
   - [Cosine Similarity (Cos)](#cosine-similarity-cos)
   - [Spectral Entropy Similarity (Ent)](#spectral-entropy-similarity-ent)
   - [Dice Coefficient (Dice)](#dice-coefficient-dice)
   - [Jaccard Index (Jacc)](#jaccard-index-jacc)
   - [Shared Peak Ratio (ShR)](#shared-peak-ratio-shr)
3. [Comparison Summary Table](#comparison-summary-table)
4. [Context Columns](#context-columns)
5. [17 Ion Match Statistics](#17-ion-match-statistics)
6. [Interpretation Guide](#interpretation-guide)

---

## Overview

After annotating both query and reference spectra with theoretical fragment ions
(b, y, etc.) within a configured m/z tolerance, MSpecViz computes two categories
of similarity metrics:

**Annotated-peak metrics** consider only peaks matched to a theoretical ion.
These measure how well the identified fragments agree between two spectra.

**All-peak metrics** consider every peak in the spectrum, whether annotated or
not. These capture overall spectral shape including noise and unannotated
fragments.

Using both categories together distinguishes spectra that share the same
fragment identities (high annotated metrics) from spectra that are globally
similar in shape (high all-peak metrics).

### Terminology

Given two annotated spectra Q (query) and R (reference):

- **A_Q** = set of annotation labels in Q where `is_annotated == True`
- **A_R** = set of annotation labels in R where `is_annotated == True`
- **S** = A_Q intersect A_R (shared annotations present in both spectra)
- For each shared label _l_ in S: **q_l** and **r_l** are the intensities
  of ion _l_ in Q and R respectively

---

## 7 Comparison Metrics

### R-squared Matching (R2m)

**Category:** Annotated-peak

**What it measures:** How well a linear model predicts reference intensities
from query intensities, using only the ions that appear in both spectra.

**Method:**
1. Find all shared annotation labels S = A_Q intersect A_R.
2. Build intensity vectors q = (q_l1, q_l2, ..., q_ln) and
   r = (r_l1, r_l2, ..., r_ln).
3. Fit a linear regression: r = a * q + b.
4. Return the Pearson R-squared (coefficient of determination).

**Formula:**

```
R2 = 1 - sum((r_i - r_hat_i)^2) / sum((r_i - r_bar)^2)
```

where r_hat_i are the predicted values from the fitted model and r_bar is the
mean of the reference intensities.

**Range:** (-inf, 1]. A value of 1.0 indicates a perfect linear relationship.
Negative values indicate a model worse than a horizontal mean line.

**Minimum data:** Requires at least 2 shared ions. Returns NaN otherwise.

**When to use:** To assess whether intensity ratios are preserved between
spectra. Tolerant of different absolute intensity scales because the regression
intercept absorbs any offset. Best for evaluating quantitative agreement among
shared fragment ions.

---

### R-squared Binned (R2b)

**Category:** All-peak

**What it measures:** How well the overall intensity profiles of two spectra
correlate when discretized into fixed-width m/z bins.

**Method:**
1. Determine the combined m/z range across both spectra.
2. Divide this range into 1000 equal-width bins.
3. For each spectrum, sum peak intensities falling into each bin.
4. Keep only bins where both spectra have signal (intersection).
5. Fit a linear regression on the bin-intensity pairs.
6. Return R-squared.

**Range:** (-inf, 1]. Typically [0, 1] for similar spectra.

**Minimum data:** Requires at least 2 shared bins. Returns NaN otherwise.

**When to use:** To compare overall spectral shape without requiring
annotation. Useful when annotation is incomplete or unreliable, or when you
want to capture unannotated features. Robust to small m/z shifts because
nearby peaks contribute to the same or adjacent bins.

---

### Cosine Similarity (Cos)

**Category:** Annotated-peak

**What it measures:** The angular similarity between the intensity vectors of
shared annotated ions.

**Formula:**

```
Cos = (q dot r) / (||q|| * ||r||)
    = sum(q_i * r_i) / (sqrt(sum(q_i^2)) * sqrt(sum(r_i^2)))
```

**Range:** [0, 1] for non-negative intensities (which MS spectra always have).
1.0 means the intensity vectors point in exactly the same direction.

**Minimum data:** Requires at least 1 shared ion. Returns NaN otherwise.

**When to use:** The de facto standard for spectral library searching (used by
NIST, MassBank). Scale-invariant: multiplying all intensities by a constant
does not change the result. Note that it only considers shared ions -- use
Dice or Jaccard alongside it to assess how many ions actually overlap.

---

### Spectral Entropy Similarity (Ent)

**Category:** All-peak

**What it measures:** Information-theoretic similarity based on Shannon entropy.
Measures how much information is gained or lost when two spectra are merged.

**Method:**
1. Normalize both spectra so intensities sum to 1 (probability distributions).
2. Merge spectra by matching peaks within 0.02 Da tolerance and averaging
   matched intensities; unmatched peaks kept at half intensity.
3. Compute Shannon entropy of Q, R, and the merged spectrum M.

**Formula:**

```
Ent = 1 - (2 * H(M) - H(Q) - H(R)) / ln(4)
```

where H(X) = -sum(p_i * ln(p_i)) is the Shannon entropy of spectrum X after
normalizing intensities to a probability distribution.

**Range:** [0, 1]. 1.0 when spectra are identical; 0.0 when they share no peaks.

**Minimum data:** Requires at least 1 peak in each spectrum. Returns NaN otherwise.

**When to use:** Shown to outperform cosine similarity for small-molecule
compound identification in benchmarks (Li et al., Nature Methods, 2021).
Considers all peaks and handles spectra of different sizes naturally.
Implemented via the optimized `ms_entropy` C library. Requires the `ms_entropy`
package.

---

### Dice Coefficient (Dice)

**Category:** Annotated-peak

**What it measures:** The fraction of annotated ions that are shared between
both spectra, weighted symmetrically.

**Formula:**

```
Dice = 2 * |S| / (|A_Q| + |A_R|)
```

**Range:** [0, 1]. 1.0 means every annotated ion in Q also appears in R and
vice versa. 0.0 means no overlap.

**When to use:** A simple, intuitive measure of annotation overlap. Symmetric:
Dice(Q, R) = Dice(R, Q). Equivalent to the harmonic mean of precision and
recall. Does not use intensity information -- a barely-visible shared peak
counts the same as the base peak. Use alongside Cosine or R2m for a complete
picture.

---

### Jaccard Index (Jacc)

**Category:** Annotated-peak

**What it measures:** Same concept as Dice but with a stricter denominator
(union instead of sum).

**Formula:**

```
Jacc = |S| / (|A_Q| + |A_R| - |S|)
```

**Range:** [0, 1]. Always less than or equal to Dice for the same data.

**Relationship to Dice:**

```
Dice = 2 * Jacc / (1 + Jacc)
```

**When to use:** More conservative than Dice -- penalizes non-shared ions more
heavily. A standard set-similarity metric. Like Dice, ignores intensity.

---

### Shared Peak Ratio (ShR)

**Category:** Annotated-peak

**What it measures:** What fraction of the smaller annotated set is shared.

**Formula:**

```
ShR = |S| / min(|A_Q|, |A_R|)
```

**Range:** [0, 1]. 1.0 means every annotated ion in the smaller spectrum has a
match in the larger one.

**When to use:** Useful when spectra have very different numbers of annotated
peaks (for example, a low-energy query with few fragments versus a high-energy
reference with many). Answers the question: "does the smaller spectrum's
information fully appear in the larger one?" Note that a query with only 1
annotated ion matching can yield ShR = 1.0, so always check the context
columns for absolute counts.

---

## Comparison Summary Table

| Property | R2m | R2b | Cos | Ent | Dice | Jacc | ShR |
|---|---|---|---|---|---|---|---|
| Category | Annotated | All-peak | Annotated | All-peak | Annotated | Annotated | Annotated |
| Uses intensity | Yes | Yes | Yes | Yes | No | No | No |
| Range | (-inf, 1] | (-inf, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] |
| Perfect match | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| No overlap | NaN | NaN | NaN | 0.0 | 0.0 | 0.0 | NaN |
| Scale-invariant | Yes | Yes | Yes | ~Yes | N/A | N/A | N/A |
| Min. data | 2 shared | 2 shared bins | 1 shared | 1 peak each | 1 ann. each | 1 ann. each | 1 ann. each |

---

## Context Columns

The comparison output includes context columns alongside the 7 metrics:

| Column | Description |
|---|---|
| Peaks (Q) | Total number of peaks in the query spectrum |
| Peaks (R) | Total number of peaks in the reference spectrum |
| Annotated (Q) | Number of peaks in Q matched to a theoretical ion |
| Annotated (R) | Number of peaks in R matched to a theoretical ion |
| Shared Ions | Number of annotation labels present in both Q and R |

These counts explain why a metric might be NaN or unexpectedly low:

- R2m = NaN means Shared Ions < 2 (not enough points for regression).
- Cos = NaN means Shared Ions = 0.
- Low Dice with high Cos means few shared ions, but those that match have
  very similar intensities.
- High Dice with low R2m means many shared ion labels, but intensity ratios
  are poorly preserved.

---

## 17 Ion Match Statistics

These supplementary metrics are computed by PyOpenMS `SpectrumAnnotator.addIonMatchStatistics()` for each query spectrum. They provide detailed fragment ion analysis beyond the pairwise comparison metrics above. In the output, all column names are prefixed with `ims_`.

### Current and Intensity Statistics

| Statistic | Column | Description | Range |
|---|---|---|---|
| N-terminal Ion Current Ratio | `ims_NTermIonCurrentRatio` | Fraction of total ion current carried by N-terminal (a, b, c) ions | [0, 1] |
| C-terminal Ion Current Ratio | `ims_CTermIonCurrentRatio` | Fraction of total ion current carried by C-terminal (x, y, z) ions | [0, 1] |
| Total Ion Current | `ims_TIC` | Sum of all peak intensities in the spectrum | [0, inf) |
| Max Ion Current | `ims_MIC` | Intensity of the most intense peak | [0, inf) |
| Annotated TIC Fraction | `ims_annotated_tic_fraction` | Fraction of total ion current explained by annotated (matched) peaks | [0, 1] |

### Fragment Ion Coverage

| Statistic | Column | Description | Range |
|---|---|---|---|
| Matched Prefix Ions Fraction | `ims_matched_prefix_ions_fraction` | Fraction of possible N-terminal fragment ions that were matched | [0, 1] |
| Matched Suffix Ions Fraction | `ims_matched_suffix_ions_fraction` | Fraction of possible C-terminal fragment ions that were matched | [0, 1] |
| Coverage | `ims_coverage` | Overall backbone cleavage coverage (fraction of peptide bonds with at least one matched fragment) | [0, 1] |
| Annotated Ion Fraction | `ims_annotated_ion_fraction` | Fraction of all peaks in the spectrum that were annotated | [0, 1] |

### Series Analysis

| Statistic | Column | Description | Range |
|---|---|---|---|
| Max Series Type | `ims_max_series_type` | Ion type of the longest consecutive series (e.g., "y") | categorical |
| Max Series Charge | `ims_max_series_charge` | Charge state of the longest consecutive series | integer |
| Max Series Size | `ims_max_series_size` | Length of the longest consecutive ion series | integer |

### Mass Accuracy

| Statistic | Column | Description | Range |
|---|---|---|---|
| MS2 Scan Index | `ims_ms2_scan_index` | Index of the MS2 scan in the original file | integer |
| Precursor m/z Error (ppm) | `ims_precursor_mz_error_ppm` | Difference between observed and theoretical precursor m/z in ppm | (-inf, inf) |
| Fragment Mass Error Mean (ppm) | `ims_fragment_mass_error_mean_ppm` | Mean mass error across all matched fragment ions in ppm | (-inf, inf) |
| Fragment Mass Error Stdev (ppm) | `ims_fragment_mass_error_stdev_ppm` | Standard deviation of fragment mass errors in ppm | [0, inf) |
| Median Fragment Error (ppm) | `ims_median_fragment_error_ppm` | Median mass error across all matched fragment ions in ppm | (-inf, inf) |

### Interpreting IMS Statistics

**Spectrum quality indicators:**
- High `ims_coverage` (> 0.7) indicates good fragmentation with most backbone
  bonds represented.
- High `ims_annotated_tic_fraction` (> 0.5) means most of the signal is
  explained by known fragment ions rather than noise or unknown fragments.
- Low `ims_annotated_ion_fraction` with high `ims_annotated_tic_fraction`
  indicates many low-intensity noise peaks but strong signal from real fragments.

**Mass accuracy indicators:**
- `ims_fragment_mass_error_mean_ppm` near 0 indicates well-calibrated data.
- Large `ims_fragment_mass_error_stdev_ppm` suggests inconsistent mass accuracy
  or incorrect matches.
- `ims_precursor_mz_error_ppm` helps verify precursor assignment quality.

**Fragmentation pattern indicators:**
- `ims_NTermIonCurrentRatio` vs `ims_CTermIonCurrentRatio` reveals whether
  fragmentation favors N-terminal or C-terminal ions.
- `ims_max_series_size` indicates the longest unbroken ion series, which
  correlates with confident sequence assignment.

---

## Interpretation Guide

### Choosing the Right Metric

**"Do these spectra contain the same fragment ions?"**
Use **Dice** or **Jaccard**. These are purely set-based and report overlap
regardless of intensity.

**"Are the shared fragments at similar relative intensities?"**
Use **R2m** or **Cosine**. These measure intensity agreement among shared ions.

**"Do the spectra have the same overall shape?"**
Use **R2b** or **Entropy**. These consider all peaks and are useful when
annotation is incomplete.

**"Is the smaller spectrum fully contained in the larger one?"**
Use **Shared Peak Ratio (ShR)**. This answers whether the information in the
less-rich spectrum is fully captured.

### Common Patterns

| R2m | R2b | Cos | Ent | Dice | Interpretation |
|---|---|---|---|---|---|
| High | High | High | High | High | Excellent match across all measures |
| High | Low | High | Low | High | Annotated ions agree, but background/noise differs |
| Low | High | Low | High | Low | Overall shape is similar, but specific ion assignments differ |
| NaN | High | NaN | High | 0.0 | No shared annotations but spectra look alike -- possible annotation issue |
| High | High | High | High | Low | Few shared ions, but those that match are in good agreement -- likely incomplete fragmentation |

---

## References

1. Stein, S.E. and Scott, D.R. (1994). "Optimization and Testing of Mass
   Spectral Library Search Algorithms for Compound Identification." Journal
   of the American Society for Mass Spectrometry, 5(9), 859-866.

2. Li, Y. et al. (2021). "Spectral entropy outperforms MS/MS dot product
   similarity for small-molecule compound identification." Nature Methods,
   18, 1524-1531. DOI: 10.1038/s41592-021-01331-z
