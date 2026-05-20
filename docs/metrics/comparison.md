# 7 Comparison Metrics

Quantitative similarity measures computed for every query-reference pair.
Four intensity-aware (R2m, R2b, Cos, Ent), three set-based (Dice, Jacc,
ShR). Implementation:

- Intensity helpers (matched-array R², binned R², matched-array cosine)
  live in `core/annotation/intensity_metrics.py`.
- Set-overlap metrics, the entropy wrapper, and the public aggregator
  `compute_all_metrics()` live in `core/annotation/comparison_metrics.py`.

Both modules are re-exported from `core.annotation`.

For a decision table ("which metric should I use?"), see
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

**What it measures:** How tightly the shared-annotation intensities of
query and reference fall on a straight line.

**Method:**
1. Group peaks by annotation label inside each spectrum and keep the
   brightest peak per label (collapsing isotopes / charge variants).
2. Find shared labels S = A_Q ∩ A_R.
3. Build intensity vectors q, r over those shared labels.
4. Compute the squared Pearson correlation coefficient via
   `scipy.stats.linregress(q, r).rvalue ** 2`.

**Formula:**

```
R² = corr(q, r)²
   = (Σ(q_i − q̄)(r_i − r̄))² / (Σ(q_i − q̄)² · Σ(r_i − r̄)²)
```

**Range:** [0, 1]. 1.0 = perfectly collinear shared intensities. Because
we use the squared correlation (not the regression coefficient of
determination), the value never goes negative even when slope < 0.

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
1. Determine the combined m/z range across both spectra (floor of joint
   min to ceil of joint max).
2. Divide this range into 1 000 equal-width bins.
3. For each spectrum, sum peak intensities into each bin via
   `numpy.bincount` over `numpy.digitize`.
4. Keep only bins where both spectra have signal (intersection).
5. Return the squared Pearson correlation coefficient on the surviving
   bin-intensity pairs.

**Range:** [0, 1] (squared correlation, never negative).

**Minimum data:** Either spectrum empty or degenerate m/z range → NaN.
≥ 2 shared bins required; returns NaN otherwise.

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
direction. Implemented as `1 − scipy.spatial.distance.cosine(q, r)`.

**Minimum data:** ≥ 1 shared ion **and** non-zero norm on both vectors.
Returns NaN when no labels are shared or when one vector is all-zero
(which can happen if the brightest peak per shared label is 0 after
normalization — rare in practice).

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

**Method:** Delegated to `ms_entropy.calculate_entropy_similarity()` with
`ms2_tolerance_in_da` set to the configured `mz_tolerance` (default
0.02 Da, from `DEFAULTS.mz_tolerance`). The library internally:

1. Normalizes both spectra to probability distributions.
2. Matches peaks within tolerance and merges into a combined spectrum M.
3. Returns the Tsallis/Shannon-style similarity below.

**Formula:**

```
Ent = 1 − (2·H(M) − H(Q) − H(R)) / ln(4)
```

where H(X) = −Σ pᵢ·ln(pᵢ) is the Shannon entropy of X after normalizing
intensities to a probability distribution.

**Range:** [0, 1]. 1.0 when spectra are identical; near 0 when they
share no peaks.

**Minimum data:** Both spectra non-empty. Returns NaN if either is empty,
if `ms_entropy` isn't installed (optional dependency), or if the library
raises (`ValueError`, `FloatingPointError`) — those errors are downgraded
to NaN so callers can tabulate.

**When to use:** Shown to outperform cosine similarity for small-molecule
compound identification in benchmarks (Li et al., Nature Methods, 2021).
Considers all peaks and handles spectra of different sizes naturally.
The tolerance is the same `mz_tolerance` used for fragment annotation —
tweak both together for consistent matching.

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

**Minimum data:** `|A_Q| + |A_R| > 0` (i.e. at least one of the two
spectra has any annotated ions). Returns NaN only when both spectra are
fully unannotated.

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

**Minimum data:** `|A_Q| + |A_R| − |S| > 0`. Returns NaN only when both
spectra are fully unannotated (in which case S = 0 too).

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

**Minimum data:** `min(|A_Q|, |A_R|) > 0` — at least one annotated ion
in each spectrum. Returns NaN when either spectrum has zero annotations.

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
| Range | [0, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] | [0, 1] |
| Perfect match | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| No overlap | NaN | NaN | NaN | ~0.0 | 0.0 | 0.0 | 0.0 |
| Scale-invariant | Yes | Yes | Yes | ~Yes | N/A | N/A | N/A |
| Min. data | 2 shared | 2 shared bins | 1 shared, non-zero norm | both non-empty | any annotation | any annotation | ≥1 ann. each |

---

## Context columns

The comparison output includes context columns alongside the 7 metrics:

| Column key | Description |
|---|---|
| `query_peak_count` | Total number of peaks in the query spectrum |
| `reference_peak_count` | Total number of peaks in the reference spectrum |
| `query_annotation_count` | **Unique** annotation labels in Q (duplicates within Q collapse to one label) |
| `reference_annotation_count` | **Unique** annotation labels in R |
| `shared` | Number of annotation labels present in both Q and R (\|S\|) |

These counts explain why a metric might be NaN or unexpectedly low:

- R2m = NaN means Shared Ions < 2 (not enough points for regression).
- Cos = NaN means Shared Ions = 0.
- Low Dice with high Cos means few shared ions, but those that match have
  very similar intensities.
- High Dice with low R2m means many shared ion labels, but intensity
  ratios are poorly preserved.
