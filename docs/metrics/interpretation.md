# Choosing and Interpreting Metrics

Once you have all 7 comparison metrics + 17 ion match statistics, which
one answers your question? This page is the decision guide.

For formulas see [comparison.md](comparison.md) and
[ion-match-stats.md](ion-match-stats.md).

---

## Choosing the right metric

**"Do these spectra contain the same fragment ions?"**
Use **Dice** or **Jaccard**. These are purely set-based and report
overlap regardless of intensity.

**"Are the shared fragments at similar relative intensities?"**
Use **R2m** or **Cosine**. These measure intensity agreement among
shared ions.

**"Do the spectra have the same overall shape?"**
Use **R2b** or **Entropy**. These consider all peaks and are useful when
annotation is incomplete.

**"Is the smaller spectrum fully contained in the larger one?"**
Use **Shared Peak Ratio (ShR)**. This answers whether the information in
the less-rich spectrum is fully captured.

**"Is the query a well-fragmented, high-quality spectrum?"**
Look at `ims_coverage` and `ims_annotated_tic_fraction`. Together they
tell you both how many backbone bonds broke and how much of the intensity
is explained.

**"Is the instrument calibrated correctly for this run?"**
Look at `ims_fragment_mass_error_mean_ppm` and
`ims_precursor_mz_error_ppm`. Near-zero with low standard deviation =
good calibration.

---

## Common patterns

| R2m | R2b | Cos | Ent | Dice | Interpretation |
|---|---|---|---|---|---|
| High | High | High | High | High | Excellent match across all measures |
| High | Low | High | Low | High | Annotated ions agree, but background/noise differs |
| Low | High | Low | High | Low | Overall shape is similar, but specific ion assignments differ |
| NaN | High | NaN | High | 0.0 | No shared annotations but spectra look alike — possible annotation issue |
| High | High | High | High | Low | Few shared ions, but those that match are in good agreement — likely incomplete fragmentation |

---

## Cross-checks

Always sanity-check with the context columns (Peaks (Q), Peaks (R),
Shared Ions):

- **ShR = 1.0 but only 1 shared ion.** Technically perfect, practically
  meaningless. Require a minimum of (say) 5 shared ions before trusting
  ShR.
- **R2m near 0 but Cos near 1.** Cos is scale-invariant; R² is not.
  Check whether the reference has been globally rescaled (`intensity_mode`
  mismatch).
- **Entropy low but Dice high.** A handful of high-intensity shared ions
  dominates entropy while the sets agree. Peak count may be small.

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
- [../architecture/core.md#mspecvizannotation](../architecture/core.md#mspecvizannotation)
  — where each metric is implemented in the codebase.
