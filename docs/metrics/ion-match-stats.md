# 17 Ion Match Statistics

Supplementary per-spectrum metrics from PyOpenMS
`SpectrumAnnotator.addIonMatchStatistics()`, computed for each query
spectrum during plot generation (`mspecviz/plotting/worker_execution.py`).

In the exported CSV/XLSX and in the Comparison Metrics table, every
column below is prefixed with `ims_`. These statistics describe the
*query* spectrum (not the pair), so they're orthogonal to the 7
comparison metrics in [comparison.md](comparison.md).

---

## Current and intensity statistics

| Statistic | Column | Description | Range |
|---|---|---|---|
| N-terminal Ion Current Ratio | `ims_NTermIonCurrentRatio` | Fraction of total ion current carried by N-terminal (a, b, c) ions | [0, 1] |
| C-terminal Ion Current Ratio | `ims_CTermIonCurrentRatio` | Fraction of total ion current carried by C-terminal (x, y, z) ions | [0, 1] |
| Total Ion Current | `ims_TIC` | Sum of all peak intensities in the spectrum | [0, ∞) |
| Max Ion Current | `ims_MIC` | Intensity of the most intense peak | [0, ∞) |
| Annotated TIC Fraction | `ims_annotated_tic_fraction` | Fraction of total ion current explained by annotated (matched) peaks | [0, 1] |

---

## Fragment ion coverage

| Statistic | Column | Description | Range |
|---|---|---|---|
| Matched Prefix Ions Fraction | `ims_matched_prefix_ions_fraction` | Fraction of possible N-terminal fragment ions that were matched | [0, 1] |
| Matched Suffix Ions Fraction | `ims_matched_suffix_ions_fraction` | Fraction of possible C-terminal fragment ions that were matched | [0, 1] |
| Coverage | `ims_coverage` | Overall backbone cleavage coverage (fraction of peptide bonds with at least one matched fragment) | [0, 1] |
| Annotated Ion Fraction | `ims_annotated_ion_fraction` | Fraction of all peaks in the spectrum that were annotated | [0, 1] |

---

## Series analysis

| Statistic | Column | Description | Range |
|---|---|---|---|
| Max Series Type | `ims_max_series_type` | Ion type of the longest consecutive series (e.g., "y") | categorical |
| Max Series Charge | `ims_max_series_charge` | Charge state of the longest consecutive series | integer |
| Max Series Size | `ims_max_series_size` | Length of the longest consecutive ion series | integer |

---

## Mass accuracy

| Statistic | Column | Description | Range |
|---|---|---|---|
| MS2 Scan Index | `ims_ms2_scan_index` | Index of the MS2 scan in the original file | integer |
| Precursor m/z Error (ppm) | `ims_precursor_mz_error_ppm` | Difference between observed and theoretical precursor m/z in ppm | (−∞, ∞) |
| Fragment Mass Error Mean (ppm) | `ims_fragment_mass_error_mean_ppm` | Mean mass error across all matched fragment ions in ppm | (−∞, ∞) |
| Fragment Mass Error Stdev (ppm) | `ims_fragment_mass_error_stdev_ppm` | Standard deviation of fragment mass errors in ppm | [0, ∞) |
| Median Fragment Error (ppm) | `ims_median_fragment_error_ppm` | Median mass error across all matched fragment ions in ppm | (−∞, ∞) |

---

## Interpreting IMS statistics

**Spectrum quality indicators:**
- High `ims_coverage` (> 0.7) indicates good fragmentation with most
  backbone bonds represented.
- High `ims_annotated_tic_fraction` (> 0.5) means most of the signal is
  explained by known fragment ions rather than noise or unknown
  fragments.
- Low `ims_annotated_ion_fraction` with high
  `ims_annotated_tic_fraction` indicates many low-intensity noise peaks
  but strong signal from real fragments.

**Mass accuracy indicators:**
- `ims_fragment_mass_error_mean_ppm` near 0 indicates well-calibrated
  data.
- Large `ims_fragment_mass_error_stdev_ppm` suggests inconsistent mass
  accuracy or incorrect matches.
- `ims_precursor_mz_error_ppm` helps verify precursor assignment
  quality.

**Fragmentation pattern indicators:**
- `ims_NTermIonCurrentRatio` vs `ims_CTermIonCurrentRatio` reveals
  whether fragmentation favors N-terminal or C-terminal ions.
- `ims_max_series_size` indicates the longest unbroken ion series, which
  correlates with confident sequence assignment.
