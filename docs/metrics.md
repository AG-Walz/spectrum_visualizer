# Comparison Metrics

Seven similarity measures computed for every query-reference pair in
`core/comparison_metrics.py` (`compute_all_metrics`). Four are intensity-aware
(R2m, R2b, Cos, Ent), three are set-based (Dice, Jacc, ShR). Annotated-peak metrics
use only peaks matched to a theoretical ion; all-peak metrics use every peak.
Using both kinds together separates *shared fragment identities* from *overall
spectral shape*.

**Terminology.** For annotated spectra Q (query) and R (reference): **A_Q**, **A_R**
are the sets of annotation labels in each; **S = A_Q ∩ A_R** the shared labels; and
**q_l**, **r_l** the intensities of a shared label *l*.

---

## The seven metrics

**R2m — R-squared matching** *(annotated)*. Squared Pearson correlation of the
shared-label intensities (brightest peak per label), via `scipy.stats.linregress`
(squared `rvalue`). Range [0, 1] (squared → never negative).
Needs ≥2 shared ions, else NaN. Best for: are intensity *ratios* preserved? Tolerant
of a global rescale (the intercept absorbs offsets).

```
R² = corr(q, r)² = (Σ(qᵢ−q̄)(rᵢ−r̄))² / (Σ(qᵢ−q̄)²·Σ(rᵢ−r̄)²)
```

**R2b — R-squared binned** *(all-peak)*. Squared Pearson correlation over a
fixed-width m/z grid spanning the joint range (`n_bins`, default 1000), summing
intensities per bin (`numpy.histogram`) and keeping bins where both have signal.
Range [0, 1]. Needs ≥2 shared bins, else NaN. Best for: overall shape without
annotation; robust to small m/z shifts (nearby peaks share a bin).

**Cos — cosine similarity** *(annotated)*. Angular similarity of shared-ion
intensity vectors, via `1 − scipy.spatial.distance.cosine`. Range [0, 1] for
non-negative MS intensities. Needs ≥1 shared ion and non-zero norm on both, else NaN.
The de facto library-search standard (NIST, MassBank); scale-invariant. Pair with
Dice/Jaccard to see how *many* ions overlap.

```
Cos = (q·r) / (‖q‖·‖r‖) = Σ(qᵢ·rᵢ) / (√Σqᵢ²·√Σrᵢ²)
```

**Ent — spectral entropy similarity** *(all-peak)*. Information-theoretic; delegated
to `ms_entropy.calculate_entropy_similarity()` with `ms2_tolerance_in_da` =
`mz_tolerance`. Range [0, 1]. NaN if either spectrum is empty, `ms_entropy` is
missing (optional dep), or the library raises. Outperforms cosine for small-molecule
ID (Li et al. 2021); shares the annotation tolerance, so tune both together.

```
Ent = 1 − (2·H(M) − H(Q) − H(R)) / ln(4),   H(X) = −Σ pᵢ·ln(pᵢ)
```

(M = the tolerance-merged spectrum; H is Shannon entropy of the normalised intensities.)

**Dice coefficient** *(annotated, set-based)*. Symmetric fraction of shared annotated
ions. Range [0, 1]; NaN only when both spectra are fully unannotated. Equals the
harmonic mean of precision and recall. Ignores intensity.

```
Dice = 2·|S| / (|A_Q| + |A_R|)
```

**Jacc — Jaccard index** *(annotated, set-based)*. Same idea, stricter denominator
(union). Always ≤ Dice; NaN only when both spectra are fully unannotated.

```
Jacc = |S| / (|A_Q| + |A_R| − |S|),    Dice = 2·Jacc / (1 + Jacc)
```

**ShR — shared peak ratio** *(annotated, set-based)*. Fraction of the *smaller*
annotated set that is shared. Range [0, 1]; NaN when either spectrum has zero
annotations. Useful when the two spectra have very different annotated-peak counts —
but a single matching ion yields ShR = 1.0, so always check the `shared` count.

```
ShR = |S| / min(|A_Q|, |A_R|)
```

| Property | R2m | R2b | Cos | Ent | Dice | Jacc | ShR |
|---|---|---|---|---|---|---|---|
| Category | Annot. | All-peak | Annot. | All-peak | Annot. | Annot. | Annot. |
| Uses intensity | Yes | Yes | Yes | Yes | No | No | No |
| No overlap | NaN | NaN | NaN | ~0.0 | 0.0 | 0.0 | 0.0 |
| Scale-invariant | Yes | Yes | Yes | ~Yes | — | — | — |
| Min. data | 2 shared | 2 bins | 1 shared, ≠0 norm | both non-empty | any annot. | any annot. | ≥1 each |

All range [0, 1] with a perfect match at 1.0.

---

## Context columns

`compute_all_metrics` also emits `query_peak_count`, `reference_peak_count`,
`query_annotation_count`, `reference_annotation_count` (unique labels), and `shared`
(|S|) — these ship as the context-count columns of the comparison summary
([export.md → data files](export.md#data-files)). They explain a NaN or surprising value:

- R2m = NaN → shared ions < 2; Cos = NaN → shared ions = 0.
- Low Dice + high Cos → few shared ions, but those match in intensity.
- High Dice + low R2m → many shared labels, but intensity ratios diverge.

---

## Choosing a metric

| Question | Use |
|---|---|
| Same fragment ions present? | Dice / Jaccard (set overlap, intensity-blind) |
| Shared fragments at similar intensities? | R2m / Cosine |
| Same overall shape? | R2b / Entropy (all-peak; annotation-free) |
| Smaller spectrum contained in larger? | ShR (and check `shared`) |
| Enough fragmentation to trust the match? | `query_annotation_count`, `shared` |

### Common patterns

| R2m | R2b | Cos | Ent | Dice | Interpretation |
|---|---|---|---|---|---|
| High | High | High | High | High | Excellent match across the board |
| High | Low | High | Low | High | Annotated ions agree, background/noise differs |
| Low | High | Low | High | Low | Similar shape, different ion assignments |
| NaN | High | NaN | High | 0.0 | No shared annotations but spectra look alike — likely an annotation issue (wrong ion set / tolerance too tight) |
| High | High | High | High | Low | Few shared ions but in good agreement — incomplete fragmentation |
| 0.0 | High | High | High | High | Shared intensities are flat (no variance) → R2m degenerate; check for saturation/identical normalisation |

### Cross-checks

- **ShR = 1.0 on 1 shared ion** — technically perfect, practically meaningless;
  require several shared ions first.
- **R2m ≈ 0 but Cos ≈ 1** — Cos is scale-invariant, R2m uses raw intensities; check
  for an `intensity_mode` mismatch (`Relative` vs `Raw`).
- **All metrics NaN** — usually one spectrum is empty post-filtering (every peak fell
  below tolerance, or the loader returned nothing). Inspect the mirror plot (see
  [export.md → plots](export.md#plots)).

---

## References

1. Stein, S.E. & Scott, D.R. (1994). *J. Am. Soc. Mass Spectrom.* 5(9), 859–866.
2. Li, Y. et al. (2021). "Spectral entropy outperforms MS/MS dot product similarity
   for small-molecule compound identification." *Nature Methods* 18, 1524–1531.
   DOI: 10.1038/s41592-021-01331-z
