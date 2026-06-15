# Choosing and Interpreting Metrics

Once you have all 7 comparison metrics, which one answers your question?
This page is the decision guide.

For formulas and edge cases see [comparison.md](comparison.md).

---

## Choosing the right metric

**"Do these spectra contain the same fragment ions?"**
Use **Dice** or **Jaccard**. These are purely set-based and report
overlap regardless of intensity.

**"Are the shared fragments at similar relative intensities?"**
Use **R2m** or **Cosine**. These measure intensity agreement among
shared ions. R2m is robust to a global rescale (the regression fits any
linear relationship); Cosine is scale-invariant (multiplying one side by
a constant doesn't change it).

**"Do the spectra have the same overall shape?"**
Use **R2b** or **Entropy**. These consider all peaks and are useful when
annotation is incomplete. R2b is annotation-free and robust to small m/z
shifts (peaks fall into adjacent bins); Entropy is the
information-theoretic counterpart and the published recommendation for
small-molecule library searching.

**"Is the smaller spectrum fully contained in the larger one?"**
Use **Shared Peak Ratio (ShR)**. This answers whether the information in
the less-rich spectrum is fully captured. Pair it with the `shared`
context column — ShR = 1.0 with only 1 shared ion is technically
perfect, practically meaningless.

**"Is fragmentation extensive enough to trust the match?"**
Look at the `query_annotation_count` and `shared` context columns
(serialised next to the metrics). A high Cos with `shared = 2` is much
weaker evidence than the same Cos with `shared = 20`.

---

## Common patterns

| R2m | R2b | Cos | Ent | Dice | Interpretation |
|---|---|---|---|---|---|
| High | High | High | High | High | Excellent match across all measures |
| High | Low | High | Low | High | Annotated ions agree, but background/noise differs |
| Low | High | Low | High | Low | Overall shape is similar, but specific ion assignments differ |
| NaN | High | NaN | High | 0.0 | No shared annotations but spectra look alike — possible annotation issue (wrong ion-type set, or m/z tolerance too tight) |
| High | High | High | High | Low | Few shared ions, but those that match are in good agreement — likely incomplete fragmentation |
| 0.0 | High | High | High | High | Shared annotations exist but their intensities are flat (no variance) — R2m is degenerate. Check if both spectra are saturated or normalised to the same height |

---

## Cross-checks

Always sanity-check with the context columns
(`query_peak_count`, `reference_peak_count`, `query_annotation_count`,
`reference_annotation_count`, `shared`):

- **ShR = 1.0 but only 1 shared ion.** Technically perfect, practically
  meaningless. Require a minimum of (say) 5 shared ions before trusting
  ShR.
- **R2m near 0 but Cos near 1.** Cos is scale-invariant; R² is computed
  on the raw shared intensities. Check whether the reference has been
  globally rescaled (`intensity_mode` mismatch — `Relative` vs `Raw`).
- **Entropy low but Dice high.** A handful of high-intensity shared ions
  dominates entropy while the sets agree. Peak count may be small.
- **All metrics NaN.** Usually means one of the spectra is empty
  post-filtering (every peak fell below the annotation tolerance, or the
  spectrum loader returned no peaks). Inspect the mirror plot directly.

---

## References

1. Stein, S.E. and Scott, D.R. (1994). "Optimization and Testing of Mass
   Spectral Library Search Algorithms for Compound Identification."
   *Journal of the American Society for Mass Spectrometry*, 5(9),
   859–866.
2. Li, Y. et al. (2021). "Spectral entropy outperforms MS/MS dot product
   similarity for small-molecule compound identification." *Nature
   Methods*, 18, 1524–1531. DOI: 10.1038/s41592-021-01331-z

---

## See also

- [../settings/app.md#annotation](../settings/app.md#annotation) — which
  annotation knobs (ion types, m/z tolerance, neutral losses) affect
  which metrics.
- [../export/formats.md](../export/formats.md) — the exact CSV/XLSX
  column names that carry these metrics into the export.
- [../architecture/core.md](../architecture/core.md) — the
  `core/comparison_metrics.py` tour shows where each metric is
  implemented.
