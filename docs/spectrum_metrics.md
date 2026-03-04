# Spectrum Comparison Metrics

This document describes the seven quantitative metrics computed by MSpecViz when
comparing a **query** (experimental) spectrum against a **reference** spectrum in
a mirror plot. All metrics are implemented in `src/spectrum_metrics.py`.

---

## Table of Contents

1. [Overview](#overview)
2. [Metric Categories](#metric-categories)
3. [Annotated-Peak Metrics](#annotated-peak-metrics)
   - [R² Matching (R2m)](#r²-matching-r2m)
   - [Cosine Similarity (Cos)](#cosine-similarity-cos)
   - [Dice Coefficient (Dice)](#dice-coefficient-dice)
   - [Jaccard Index (Jacc)](#jaccard-index-jacc)
   - [Shared Peak Ratio (ShR)](#shared-peak-ratio-shr)
4. [All-Peak Metrics](#all-peak-metrics)
   - [R² Binned (R2b)](#r²-binned-r2b)
   - [Spectral Entropy Similarity (Ent)](#spectral-entropy-similarity-ent)
5. [Context Columns](#context-columns)
6. [Comparison Table](#comparison-table)
7. [Interpretation Guide](#interpretation-guide)
8. [References](#references)

---

## Overview

After annotating both the query and reference spectra with theoretical fragment
ions (b, y, etc.) within a given m/z tolerance, MSpecViz computes similarity
metrics that fall into two complementary categories:

- **Annotated-peak metrics** — only consider peaks that were successfully
  matched to a theoretical ion. These measure how well the *identified*
  fragments agree between the two spectra.

- **All-peak metrics** — consider every peak in the spectrum, whether annotated
  or not. These capture the overall spectral shape including noise, unannotated
  fragments, and background.

Using both categories together lets you distinguish between spectra that share
the same fragment identities (high annotated metrics) versus spectra that are
globally similar in shape (high all-peak metrics). A spectrum pair with high
annotated metrics but low all-peak metrics would indicate that the known
fragments match but the background/noise patterns differ significantly.

---

## Metric Categories

| Metric | Abbreviation | Category | Peaks Used |
|--------|-------------|----------|------------|
| R² Matching | R2m | Annotated-peak | Shared annotated ions only |
| Cosine Similarity | Cos | Annotated-peak | Shared annotated ions only |
| Dice Coefficient | Dice | Annotated-peak | All annotated ions (shared & unshared) |
| Jaccard Index | Jacc | Annotated-peak | All annotated ions (shared & unshared) |
| Shared Peak Ratio | ShR | Annotated-peak | All annotated ions (shared & unshared) |
| R² Binned | R2b | All-peak | Every peak (binned across m/z range) |
| Spectral Entropy | Ent | All-peak | Every peak (m/z tolerance-based matching) |

---

## Annotated-Peak Metrics

These metrics operate on peaks that have been assigned a theoretical ion label
(e.g., `b3`, `y7++`). They require a prior annotation step.

### Terminology

Given two annotated spectra Q (query) and R (reference):

- **A_Q** = set of annotation labels in Q where `is_annotated == True`
- **A_R** = set of annotation labels in R where `is_annotated == True`
- **S** = A_Q ∩ A_R (shared annotations present in both spectra)
- For each shared label _l_ in S:
  - **q_l** = intensity of ion _l_ in Q
  - **r_l** = intensity of ion _l_ in R

---

### R² Matching (R2m)

**What it measures:** How well a linear model predicts reference intensities
from query intensities, using only the ions that appear in both spectra.

**Method:**
1. Find all shared annotation labels S = A_Q ∩ A_R
2. Build intensity vectors: **q** = (q_l₁, q_l₂, ..., q_lₙ) and
   **r** = (r_l₁, r_l₂, ..., r_lₙ)
3. Fit a linear regression: r = a·q + b
4. Return the coefficient of determination R²

**Equation:**

```
         Σᵢ (rᵢ - r̂ᵢ)²
R² = 1 - ───────────────
         Σᵢ (rᵢ - r̄)²
```

where r̂ᵢ are the predicted values from the fitted model and r̄ is the mean
of the reference intensities.

**Range:** (-∞, 1]. A value of 1.0 means perfect linear relationship.
Negative values indicate a model worse than a horizontal mean line.

**Minimum data:** Requires ≥ 2 shared ions. Returns NaN otherwise.

**Advantages:**
- Captures whether intensity *ratios* are preserved between spectra, not just
  presence/absence
- Tolerant of different absolute intensity scales (the intercept absorbs offset)
- Directly comparable to the specviz `spectrum_correlation()` inhouse method

**Disadvantages:**
- Ignores ions present in only one spectrum
- Sensitive to outliers when few shared ions exist
- Can be misleadingly high with only 2-3 shared ions (overfitting)

**Origin:** Adapted from the specviz `spectrum_correlation()` function.

---

### Cosine Similarity (Cos)

**What it measures:** The angular similarity between the intensity vectors of
shared annotated ions.

**Method:**
1. Find all shared annotation labels S = A_Q ∩ A_R
2. Build intensity vectors **q** and **r** (same as R²m)
3. Compute the cosine of the angle between them

**Equation:**

```
            q · r           Σᵢ qᵢ · rᵢ
Cos = ─────────────── = ─────────────────────────
      ‖q‖ · ‖r‖       √(Σᵢ qᵢ²) · √(Σᵢ rᵢ²)
```

**Range:** [0, 1] for non-negative intensities (which MS spectra always are).
1.0 means the intensity vectors point in exactly the same direction.

**Minimum data:** Requires ≥ 1 shared ion. Returns NaN if no shared ions.

**Advantages:**
- Scale-invariant: multiplying all intensities by a constant doesn't change the
  result
- The de facto standard in spectral library searching (NIST, MassBank)
- Intuitive geometric interpretation

**Disadvantages:**
- Only considers shared ions — a spectrum with 2 matching ions and 50
  non-matching ions still gets Cos=1.0 if the 2 intensities match perfectly
- Does not account for the *number* of matched peaks (that is what Dice/Jaccard
  measure)

**Origin:** Stein & Scott, 1994. "Optimization and Testing of Mass Spectral
Library Search Algorithms for Compound Identification." JASMS.

---

### Dice Coefficient (Dice)

**What it measures:** The fraction of annotated ions that are shared between
both spectra, weighted symmetrically.

**Equation:**

```
           2 · |S|
Dice = ─────────────
       |A_Q| + |A_R|
```

where |S| is the count of shared annotation labels, |A_Q| is the count of
annotated peaks in the query, and |A_R| is the count in the reference.

**Range:** [0, 1]. A value of 1.0 means every annotated ion in Q also
appears in R and vice versa. 0.0 means no overlap.

**Advantages:**
- Simple, intuitive measure of overlap
- Symmetric: Dice(Q,R) = Dice(R,Q)
- Balances between precision and recall (it is their harmonic mean)

**Disadvantages:**
- Ignores intensity information entirely — a barely-visible shared peak counts
  the same as the base peak
- Sensitive to annotation quality (missed annotations deflate the score)

---

### Jaccard Index (Jacc)

**What it measures:** Same concept as Dice but with a stricter denominator
(union instead of sum).

**Equation:**

```
              |S|
Jacc = ─────────────────
       |A_Q| + |A_R| - |S|
```

The denominator is the total number of *unique* annotation labels across both
spectra (the union).

**Range:** [0, 1]. Always ≤ Dice for the same data. Identical spectra → 1.0.

**Relationship to Dice:**

```
         2 · Jacc
Dice = ───────────
        1 + Jacc
```

**Advantages:**
- More conservative than Dice — penalizes non-shared ions more heavily
- Standard set-similarity metric used across many fields

**Disadvantages:**
- Same as Dice: ignores intensity, sensitive to annotation quality
- More conservative can also mean it undervalues partial matches

---

### Shared Peak Ratio (ShR)

**What it measures:** What fraction of the *smaller* annotated set is shared.
This is a completeness measure.

**Equation:**

```
            |S|
ShR = ─────────────────
      min(|A_Q|, |A_R|)
```

**Range:** [0, 1]. A value of 1.0 means every annotated ion in the smaller
spectrum has a match in the larger one.

**Advantages:**
- Useful when spectra have very different numbers of annotated peaks (e.g.,
  a low-energy query with few fragments vs a high-energy reference with many)
- Answers the question: "does the smaller spectrum's information fully appear
  in the larger one?"

**Disadvantages:**
- Asymmetric in spirit — dominated by the smaller set
- A query with only 1 annotated ion matching gets ShR=1.0, which can be
  misleading

---

## All-Peak Metrics

These metrics use every peak in the spectrum — annotated, unannotated, noise,
background ions, neutral losses not in the theoretical model, etc.

### R² Binned (R2b)

**What it measures:** How well the overall intensity profiles of two spectra
correlate when discretized into fixed-width m/z bins.

**Method:**
1. Determine the combined m/z range: [min(Q∪R), max(Q∪R)]
2. Divide this range into N bins (default: 1000) of equal width
3. For each spectrum, sum peak intensities falling into each bin
4. Keep only the bins where **both** spectra have signal (intersection)
5. Fit a linear regression on these bin-intensity pairs
6. Return R²

**Equation:**

Same R² formula as R2m, but applied to binned intensity vectors instead of
individual ion intensities:

```
         Σⱼ (Bʳⱼ - B̂ʳⱼ)²
R² = 1 - ──────────────────
         Σⱼ (Bʳⱼ - B̄ʳ)²
```

where j iterates over shared bins, Bʳⱼ is the summed reference intensity in
bin j, and B̂ʳⱼ is the linear prediction from the query bin intensity.

**Range:** (-∞, 1]. Usually [0, 1] for similar spectra.

**Bin intersection:** Only bins where both spectra contribute signal are
used for the regression. This follows the specviz `new_spectrum_correlation()`
convention and avoids the regression being dominated by the many empty bins
that would otherwise drag R² toward zero.

**Advantages:**
- Does not require annotation — works on raw peak lists
- Captures overall spectral shape including unannotated features
- Robust to small m/z shifts (peaks in neighboring bins still contribute to
  the same or adjacent bins)
- Useful for comparing spectra where annotation is incomplete or unreliable

**Disadvantages:**
- Bin width affects the result: too few bins → loss of detail; too many →
  sparse signals
- Intense single peaks dominate the regression
- The intersection-only approach can overestimate similarity if both spectra
  only share a few bins
- Does not distinguish between biologically meaningful and noise peaks

**Origin:** Adapted from the specviz `new_spectrum_correlation()` function.
Binning strategy references: Borgelt et al., 2016 (DOI: 10.1007/978-1-4939-6740-7_7).

---

### Spectral Entropy Similarity (Ent)

**What it measures:** Information-theoretic similarity based on Shannon
entropy. Measures how much information is gained (or lost) when two spectra
are merged.

**Method (conceptual):**
1. Normalize both spectra so intensities sum to 1 (probability distributions)
2. Merge spectra by matching peaks within an m/z tolerance (0.02 Da) and
   averaging matched intensities, keeping unmatched peaks at half intensity
3. Compute Shannon entropy of Q, R, and the merged spectrum M
4. Return the similarity score

**Equation:**

```
                2 · H(M) - H(Q) - H(R)
Ent = 1 − ─────────────────────────────
                     ln(4)
```

where:

```
H(X) = − Σᵢ pᵢ · ln(pᵢ)
```

is the Shannon entropy of spectrum X after normalizing intensities to a
probability distribution (pᵢ = Iᵢ / ΣI).

The denominator ln(4) normalizes the divergence so that Ent ∈ [0, 1]:
- **Ent = 1.0** when the spectra are identical
- **Ent = 0.0** when they share no peaks at all

**Advantages:**
- Information-theoretically grounded
- Considers all peaks, not just annotated ones
- Handles spectra of different sizes naturally
- Shown to outperform cosine similarity for spectral library matching in
  benchmarks (Li et al., 2021)
- Implemented via the optimized `ms_entropy` C library (fast)

**Disadvantages:**
- Requires the `ms_entropy` package (pip dependency)
- The m/z tolerance used for merging (0.02 Da) is fixed and may not be
  optimal for all instruments
- Less intuitive than cosine or R² — harder to explain to collaborators
- Sensitive to the intensity distribution: spectra with many low-intensity
  peaks will have higher entropy regardless of biological meaning

**Origin:** Li et al., 2021. "Spectral entropy outperforms MS/MS dot product
similarity for small-molecule compound identification." Nature Methods, 18,
1524-1531.

---

## Context Columns

In addition to the 7 metrics, the output table includes context columns that
help interpret the metric values:

| Column | Description |
|--------|-------------|
| Peaks (Q) | Total number of peaks in the query spectrum |
| Peaks (R) | Total number of peaks in the reference spectrum |
| Annotated (Q) | Number of peaks in Q matched to a theoretical ion |
| Annotated (R) | Number of peaks in R matched to a theoretical ion |
| Shared Ions | Number of annotation labels present in both Q and R |

These counts explain *why* a metric might be NaN or unexpectedly low:
- **R2m = NaN** → Shared Ions < 2 (not enough points for regression)
- **Cos = NaN** → Shared Ions = 0
- **Dice = 0.0** → Shared Ions = 0, but both spectra have annotations
- **Low Dice, high Cos** → Few shared ions, but those that match have very
  similar intensities
- **High Dice, low R2m** → Many shared ion labels, but intensity ratios are
  poorly preserved

---

## Comparison Table

| Property | R2m | R2b | Cos | Ent | Dice | Jacc | ShR |
|----------|-----|-----|-----|-----|------|------|-----|
| **Category** | Annotated | All-peak | Annotated | All-peak | Annotated | Annotated | Annotated |
| **Uses intensity** | Yes | Yes | Yes | Yes | No | No | No |
| **Range** | (-∞,1] | (-∞,1] | [0,1] | [0,1] | [0,1] | [0,1] | [0,1] |
| **Perfect match** | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| **No overlap** | NaN | NaN | NaN | 0.0 | 0.0 | 0.0 | NaN |
| **Scale-invariant** | Yes* | Yes* | Yes | ~Yes | N/A | N/A | N/A |
| **Min. data** | 2 shared | 2 shared bins | 1 shared | 1 peak each | 1 ann. each | 1 ann. each | 1 ann. each |
| **Sensitive to noise** | Low | Medium | Low | High | Low | Low | Low |
| **Sensitive to missing annotations** | High | None | High | None | High | High | High |

*R² is invariant to affine transformations (scale + offset) due to the linear
regression fit.

---

## Interpretation Guide

### When to use which metric

**"Do these spectra contain the same fragment ions?"**
→ Look at **Dice** or **Jaccard**. These are purely set-based and tell you
about overlap regardless of intensity.

**"Are the shared fragments at similar relative intensities?"**
→ Look at **R2m** or **Cosine**. These measure intensity agreement among
shared ions.

**"Do the spectra have the same overall shape?"**
→ Look at **R2b** or **Entropy**. These consider all peaks and are useful
when annotation is incomplete.

**"Is the smaller spectrum fully contained in the larger one?"**
→ Look at **SharedRatio (ShR)**. This answers whether the information in the
less-rich spectrum is fully captured.

### Common patterns

| R2m | R2b | Cos | Ent | Dice | Interpretation |
|-----|-----|-----|-----|------|---------------|
| High | High | High | High | High | Excellent match across all measures |
| High | Low | High | Low | High | Annotated ions agree, but background/noise differ |
| Low | High | Low | High | Low | Overall shape is similar, but specific ion assignments differ |
| N/A | High | N/A | High | 0.0 | No shared annotations but spectra look alike — possible annotation issue |
| High | High | High | High | Low | Few shared ions, but they match well — incomplete fragmentation |

---

## References

1. Stein, S.E. & Scott, D.R. (1994). "Optimization and Testing of Mass
   Spectral Library Search Algorithms for Compound Identification." *Journal
   of the American Society for Mass Spectrometry*, 5(9), 859-866.

2. Li, Y. et al. (2021). "Spectral entropy outperforms MS/MS dot product
   similarity for small-molecule compound identification." *Nature Methods*,
   18, 1524-1531. DOI: 10.1038/s41592-021-01331-z

3. Borgelt, C. et al. (2016). "Statistical Analysis of Proteomics, Metabolomics,
   and Lipidomics Data Using Mass Spectrometry." In: *Frontiers in Data Science*.
   Springer. DOI: 10.1007/978-1-4939-6740-7_7

4. MSpecViz internal functions: `spectrum_correlation()` and
   `new_spectrum_correlation()` in `src/spectrum_metrics.py`
