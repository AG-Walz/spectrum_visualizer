# Annotation Settings

This document explains every annotation and visualization setting in MSpecViz, what it does, and how to choose the right values for your data.

---

## Ion Types

Fragment ions are produced when a peptide bond breaks during collision-induced dissociation (CID) or similar fragmentation. Different bond positions produce different ion series.

### Available Ion Series

| Ion | Terminal | Bond Broken | Typical Fragmentation | Color |
|-----|----------|-------------|----------------------|-------|
| **b** | N-terminal | Amide C–N bond | CID, HCD | Blue (#0066CC) |
| **y** | C-terminal | Amide C–N bond | CID, HCD | Red (#CC0000) |
| **a** | N-terminal | Cα–C bond (b minus CO) | CID (minor) | Yellow (#FFCC00) |
| **c** | N-terminal | N–Cα bond + NH₃ | ETD, ECD | Green (#00CC66) |
| **x** | C-terminal | Cα–C bond + CO | ETD (rare) | Orange (#CC6600) |
| **z** | C-terminal | N–Cα bond minus NH₃ | ETD, ECD | Purple (#6600CC) |
| **precursor** | — | Intact peptide | Always present | Magenta (#8B008B) |
| **immonium** | — | Internal single-residue fragments | CID | Teal (#00CED1) |

### Recommendations

| Instrument / Method | Recommended Ion Types |
|---------------------|----------------------|
| **CID/HCD** (Orbitrap, Q-TOF, timsTOF) | b, y (default) |
| **ETD/ECD** | c, z |
| **EThcD** (hybrid) | b, y, c, z |
| **Comprehensive annotation** | b, y, a (adds minor a-ions) |

### What "Selected Ion Types" Means

When you select an ion type:
1. Theoretical fragments of that type are generated from the peptide sequence
2. Experimental peaks within the m/z tolerance of a theoretical fragment are **annotated** (colored and labeled)
3. Unmatched experimental peaks are shown in **grey** ("unannotated")

Ion types you do **not** select are simply not generated — their peaks remain grey even if they happen to fall at the right m/z.

### Precursor Ion Behavior

The precursor ion (intact peptide at various charge states) is treated specially:

- Theoretical precursor peaks are **always generated** internally (for identification purposes)
- If `precursor` is in your selected ion types → precursor peaks are **colored** (magenta) and **labeled**
- If `precursor` is not selected → precursor peaks are shown as **grey** (unannotated)
- If "Remove precursor peaks" is checked → precursor peaks are **removed entirely** from the plot

---

## m/z Tolerance

The tolerance defines how close (in Daltons) an experimental peak must be to a theoretical fragment to count as a match.

### How Matching Works

MSpecViz uses PyOpenMS's `SpectrumAlignment` algorithm:

1. Both spectra (experimental and theoretical) are sorted by m/z
2. For each theoretical peak, the algorithm searches for the closest experimental peak within the tolerance window
3. The alignment is globally optimal — it finds the assignment that maximizes matches without allowing any experimental peak to be matched twice

### Choosing a Tolerance

| Instrument Type | Typical Mass Accuracy | Recommended Tolerance |
|-----------------|----------------------|----------------------|
| timsTOF | 5–10 ppm (~0.005–0.01 Da at 1000 m/z) | 0.02 Da (default) |
| Orbitrap (high-res) | 1–5 ppm (~0.001–0.005 Da) | 0.01–0.02 Da |
| Q-TOF | 10–20 ppm | 0.02–0.05 Da |
| Triple quad / ion trap | 0.1–0.5 Da | 0.1–0.5 Da |

**The tolerance is in absolute Daltons**, not ppm. The default of 0.02 Da works well for most modern high-resolution instruments.

### Effects of Tolerance Setting

| Too tight (e.g., 0.001 Da) | Too loose (e.g., 0.5 Da) |
|----------------------------|--------------------------|
| Misses valid matches | False positive annotations |
| Low annotation coverage | Noise peaks get annotated |
| Underestimates similarity metrics | Overestimates similarity metrics |

---

## Neutral Losses

**Setting**: "Add neutral losses (H₂O, NH₃)"

When enabled, the theoretical spectrum includes peaks corresponding to fragment ions that lost a small neutral molecule:

| Loss | Mass Lost | Residues Prone to Loss |
|------|-----------|----------------------|
| H₂O | 18.011 Da | Ser, Thr, Glu, Asp |
| NH₃ | 17.027 Da | Arg, Lys, Gln, Asn |

### When to Enable

- **Enable** when you see many unannotated peaks that are systematically ~17–18 Da below annotated peaks — these are likely neutral losses
- **Disable** (default) for cleaner plots with fewer annotations. Most analyses only need the main ion series

### Effect

Enabling neutral losses roughly doubles the number of theoretical peaks. This increases annotation coverage but also increases the chance of spurious matches (especially with a loose tolerance).

---

## First Prefix Ion

**Setting**: "Add first prefix ion (a₁, b₁)"

The first prefix ion is the smallest N-terminal fragment (a single amino acid). Some search engines and theoretical generators omit it because:

- It's often below the instrument's effective scan range
- For b₁, it requires the first amide bond to break, which is sometimes disfavored

### When to Enable

- **Enable** (default) to maximize annotation coverage
- **Disable** if you're specifically comparing against a reference that also omits first prefix ions, to keep annotation sets consistent

---

## Remove Precursor Peaks

**Setting**: "Remove precursor peaks from plot"

When enabled, all peaks identified as precursor ions (at any charge state) are removed from the plot entirely. Remaining intensities are renormalized.

### When to Enable

- **Enable** when the precursor ion dominates the spectrum and obscures fragment ions. This is common with incomplete fragmentation or when the isolation window is wide.
- **Disable** (default) to see the full spectrum including the precursor.

### Behavior Details

1. Precursor peaks are identified by their annotation (annotations starting with `[M` or containing "precursor")
2. They are removed from the DataFrame before plotting
3. If intensity mode is "Relative", the base peak is recalculated after removal
4. This means previously second-tallest peaks become more visible

---

## Intensity Mode

Controls how peak heights are displayed.

### Relative (Default)

Each spectrum is independently normalized so that its tallest peak has intensity = 1.0. All other peaks are expressed as a fraction of the base peak.

**Advantages**:
- Directly comparable across spectra with different total ion currents
- Reveals the fragmentation *pattern* regardless of signal strength
- Standard for spectral library matching

**Mirror plot y-axis**: -1.0 to +1.0 (upper spectrum positive, lower negative)

### Raw

Original intensity values from the mzML file are used without normalization.

**Advantages**:
- Shows actual signal strength
- Useful for assessing spectral quality or comparing ion transmission

**Disadvantages**:
- Spectra with very different total intensities are hard to compare visually
- A faint spectrum next to an intense one will appear as barely visible lines

**Mirror plot y-axis**: Auto-scaled to the maximum intensity of either spectrum

---

## Ion Colors

Each ion type has a distinct color for easy identification in plots. The defaults are:

| Ion Type | Color | Hex |
|----------|-------|-----|
| a | Yellow | #FFCC00 |
| b | Blue | #0066CC |
| c | Green | #00CC66 |
| x | Orange | #CC6600 |
| y | Red | #CC0000 |
| z | Purple | #6600CC |
| precursor | Magenta | #8B008B |
| immonium | Teal | #00CED1 |
| unannotated | Grey | #AAAAAA |

Colors can be customized in `settings.json` under the `ion_colors` key:

```json
{
  "ion_colors": {
    "b": "#0066CC",
    "y": "#CC0000",
    "unannotated": "#AAAAAA"
  }
}
```

---

## Settings in `settings.json`

Default values for all annotation settings are defined in `settings.json`. The UI initializes from these values, and users can override them per-session.

```json
{
  "ion_types": {
    "default": ["b", "y"],
    "available": ["a", "b", "c", "x", "y", "z", "precursor", "immonium"]
  },
  "mz_tolerance": {
    "default": 0.02,
    "unit": "Da"
  },
  "annotation": {
    "add_losses": false,
    "add_first_prefix_ion": true,
    "remove_precursor_peaks": false
  }
}
```

Changes to `settings.json` take effect on the next app restart.
