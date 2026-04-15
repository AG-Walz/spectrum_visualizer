# 7 Comparison Metrics

Quantitative similarity measures computed for every query-reference pair.
Four intensity-aware, three set-based. Implemented in
`mspecviz/annotation/metrics_intensity.py` and `.../metrics_set.py`, both
re-exported from `mspecviz.annotation`.

- For the 17 ion match statistics, see [ion-match-stats.md](ion-match-stats.md).
- For a decision table ("which metric should I use?"), see
  [interpretation.md](interpretation.md).

---

## Overview

After annotating both query and reference spectra with theoretical
fragment ions (b, y, etc.) within a configured m/z tolerance, MSpecViz
computes two categories of similarity metrics:

**Annotated-peak metrics** consider only peaks matched to a theoretical
ion. These measure how well the identified fragments agree between two
spectra.

**All-peak metrics** consider every peak in the spectrum, whether
annotated or not. These capture overall spectral shape including noise
and unannotated fragments.

Using both categories together distinguishes spectra that share the same
fragment identities (high annotated metrics) from spectra that are
globally similar in shape (high all-peak metrics).

### Terminology

Given two annotated spectra Q (query) and R (reference):

- **A_Q** = set of annotation labels in Q where `is_annotated == True`
- **A_R** = set of annotation labels in R where `is_annotated == True`
- **S** = A_Q ∩ A_R (shared annotations present in both spectra)
- For each shared label *l* in S: **q_l** and **r_l** are the intensities
  of ion *l* in Q and R respectively

---

## R-squared Matching (R2m)

**Category:** Annotated-peak

**What it measures:** How well a linear model predicts reference
intensities from query intensities, using only the ions that appear in
both spectra.

**Method:**
1. Find all shared annotation labels S = A_Q ∩ A_R.
2. Build intensity vectors q = (q_l1, q_l2, …, q_ln) and
   r = (r_l1, r_l2, …, r_ln).
3. Fit a linear regression: r = a·q + b.
4. Return the Pearson R² (coefficient of determination).

**Formula:**

```
R² = 1 − Σ(r_i − r̂_i)² / Σ(r_i − r̄)²
```

**Range:** (−∞, 1]. 1.0 = perfect linear relationship. Negative values
indicate a model worse than a horizontal mean line.

**Minimum data:** ≥ 2 shared ions. Returns NaN otherwise.

**When to use:** To assess whether intensity ratios are preserved between
spectra. Tolerant of different absolute intensity scales because the
regression intercept absorbs any offset. Best for evaluating quantitative
agreement among shared fragment ions.

---

## R-squared Binned (R2b)

**Category:** All-peak

**What it measures:** How well the overall intensity profiles of two
spectra correlate when discretized into fixed-width m/z bins.

**Method:**
1. Determine the combined m/z range across both spectra.
2. Divide this range into 1 000 equal-width bins.
3. For each spectrum, sum peak intensities falling into each bin.
4. Keep only bins where both spectra have signal (intersection).
5. Fit a linear regression on the bin-intensity pairs.
6. Return R².

**Range:** (−∞, 1]. Typically [0, 1] for similar spectra.

**Minimum data:** ≥ 2 shared bins. Returns NaN otherwise.

**When to use:** To compare overall spectral shape without requiring
annotation. Useful when annotation is incomplete or unreliable, or when
you want to capture unannotated features. Robust to small m/z shifts
because nearby peaks contribute to the same or adjacent bins.

---

## Cosine Similarity (Cos)

**Category:** Annotated-peak

**What it measures:** The angular similarity between the intensity
vectors of shared annotated ions.

**Formula:**

```
Cos = (q · r) / (‖q‖ · ‖r‖)
    = Σ(q_i · r_i) / (√Σq_i² · √Σr_i²)
```

**Range:** [0, 1] for non-negative intensities (which MS spectra always
have). 1.0 means the intensity vectors point in exactly the same
direction.

**Minimum data:** ≥ 1 shared ion. Returns NaN otherwise.

**When to use:** The de facto standard for spectral library searching
(used by NIST, MassBank). Scale-invariant: multiplying all intensities by
a constant does not change the result. Note that it only considers shared
ions — use Dice or Jaccard alongside it to assess how many ions actually
overlap.

---

## Spectral Entropy Similarity (Ent)

**Category:** All-peak

**What it measures:** Information-theoretic similarity based on Shannon
entropy. Measures how much information is gained or lost when two spectra
are merged.

**Method:**
1. Normalize both spectra so intensities sum to 1 (probability
   distributions).
2. Merge spectra by matching peaks within 0.02 Da tolerance and averaging
   matched intensities; unmatched peaks kept at half intensity.
3. Compute Shannon entropy of Q, R, and the merged spectrum M.

**Formula:**

```
Ent = 1 − (2·H(M) − H(Q) − H(R)) / ln(4)
```

where H(X) = −Σ pᵢ·ln(pᵢ) is the Shannon entropy of X after normalizing
intensities to a probability distribution.

**Range:** [0, 1]. 1.0 when spectra are identical; 0.0 when they share
no peaks.

**Minimum data:** ≥ 1 peak in each spectrum. Returns NaN otherwise.

**When to use:** Shown to outperform cosine similarity for small-molecule
compound identification in benchmarks (Li et al., Nature Methods, 2021).
Considers all peaks and handles spectra of different sizes naturally.
Implemented via the optimized `ms_entropy` C library.

---

## Dice Coefficient (Dice)

**Category:** Annotated-peak

**What it measures:** The fraction of annotated ions that are shared
between both spectra, weighted symmetrically.

**Formula:**

```
Dice = 2·|S| / (|A_Q| + |A_R|)
```

**Range:** [0, 1]. 1.0 means every annotated ion in Q also appears in R
and vice versa.

**When to use:** A simple, intuitive measure of annotation overlap.
Symmetric: Dice(Q, R) = Dice(R, Q). Equivalent to the harmonic mean of
precision and recall. Does not use intensity information — a
barely-visible shared peak counts the same as the base peak. Use
alongside Cosine or R2m for a complete picture.

---

## Jaccard Index (Jacc)

**Category:** Annotated-peak

**What it measures:** Same concept as Dice but with a stricter
denominator (union instead of sum).

**Formula:**

```
Jacc = |S| / (|A_Q| + |A_R| − |S|)
```

**Range:** [0, 1]. Always ≤ Dice for the same data.

**Relationship to Dice:**

```
Dice = 2·Jacc / (1 + Jacc)
```

**When to use:** More conservative than Dice — penalizes non-shared ions
more heavily. A standard set-similarity metric. Like Dice, ignores
intensity.

---

## Shared Peak Ratio (ShR)

**Category:** Annotated-peak

**What it measures:** What fraction of the smaller annotated set is
shared.

**Formula:**

```
ShR = |S| / min(|A_Q|, |A_R|)
```

**Range:** [0, 1]. 1.0 means every annotated ion in the smaller spectrum
has a match in the larger one.

**When to use:** Useful when spectra have very different numbers of
annotated peaks (for example, a low-energy query with few fragments
versus a high-energy reference with many). Answers the question: "does
the smaller spectrum's information fully appear in the larger one?"
Note that a query with only 1 annotated ion matching can yield ShR = 1.0,
so always check the context columns for absolute counts.

---

## Comparison summary table

| Property | R2m | R2b | Cos | Ent | Dice | Jacc | ShR |
|---|---|---|---|---|---|---|---|
| Category | Annotated | All-peak | Annotated | All-peak | Annotated | Annotated | Annotated |
| Uses intensity | Yes | Yes | Yes | Yes | No | No | No |
| Range | (−∞, 1] | (−∞, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] |
| Perfect match | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| No overlap | NaN | NaN | NaN | 0.0 | 0.0 | 0.0 | NaN |
| Scale-invariant | Yes | Yes | Yes | ~Yes | N/A | N/A | N/A |
| Min. data | 2 shared | 2 shared bins | 1 shared | 1 peak each | 1 ann. each | 1 ann. each | 1 ann. each |

---

## Context columns

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
- High Dice with low R2m means many shared ion labels, but intensity
  ratios are poorly preserved.
