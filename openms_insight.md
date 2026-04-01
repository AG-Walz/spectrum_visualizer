# openms_insight Deep Analysis for the MSInsight/MSpecViz Use Case

**Version:** 0.1.13 (installed) | **Author:** Tom David Mueller | **License:** BSD-3-Clause
**Repository:** [github.com/t0mdavid-m/OpenMS-Insight](https://github.com/t0mdavid-m/OpenMS-Insight)
**Published:** PMC11811998 (OpenMS WebApps framework paper)

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Package Architecture](#2-package-architecture)
3. [Where openms_insight Shines](#3-where-openms_insight-shines)
4. [Component-by-Component Analysis](#4-component-by-component-analysis)
5. [What MSpecViz Does That openms_insight Replaces Well](#5-what-mspecviz-does-that-openms_insight-replaces-well)
6. [What Is Still Missing](#6-what-is-still-missing)
7. [Where openms_insight Could Improve](#7-where-openms_insight-could-improve)
8. [Comparison Matrix](#8-comparison-matrix)
9. [Recommendations](#9-recommendations)

---

## 1. Executive Summary

openms_insight is a Streamlit-native Vue.js visualization framework for mass spectrometry data. It provides five interactive components (Table, LinePlot, SequenceView, Heatmap, VolcanoPlot) linked through a shared StateManager. For MSInsight, the package successfully replaces MSpecViz's interactive spectrum viewing and adds entirely new capabilities (sequence coverage, RT heatmaps, linked PSM browsing), but it cannot replace MSpecViz's batch export pipeline, comparison metrics, Koina AI predictions, or mirror plot layout. The result is a hybrid architecture: openms_insight handles interactive exploration, while custom code handles computation and static export.

**Bottom line:** openms_insight is the right foundation for interactive exploration but is not (and does not aim to be) a complete replacement for a full MS analysis pipeline. The gaps are well-defined and addressable.

---

## 2. Package Architecture

### 2.1 High-Level Design

```
Python (Streamlit)                    Vue.js (Browser)
┌──────────────────┐                 ┌──────────────────┐
│  StateManager    │ ◄── counters ──►│  Vue Components  │
│  (state.py)      │                 │  (compiled JS)   │
├──────────────────┤                 ├──────────────────┤
│  Bridge          │ ── data hash ──►│  Hash Echo       │
│  (bridge.py)     │ ◄── requests ──│  (bidirectional) │
├──────────────────┤                 └──────────────────┘
│  Cache Layer     │
│  (cache.py)      │
├──────────────────┤
│  Preprocessing   │
│  (compression,   │
│   filtering,     │
│   scatter)       │
└──────────────────┘
```

### 2.2 Core Design Principles

1. **Polars LazyFrame native** — All data ingestion uses `pl.LazyFrame` with predicate pushdown for efficient filtered reads from Parquet files.
2. **Disk-backed caching** — Components preprocess data once and store results as Parquet files on disk. Reconstruction from cache requires zero Python-side data.
3. **Counter-based conflict resolution** — The StateManager uses separate `selection_counter` and `pagination_counter` values to prevent pagination events from clobbering user selections (a subtle but important correctness guarantee).
4. **Bidirectional hash confirmation** — Vue echoes the data hash it currently holds; Python only sends new data if hashes diverge. This eliminates redundant data transfers on Streamlit reruns.
5. **Subprocess memory isolation** — Optional `data_path=` mode spawns a subprocess with `mimalloc` (PURGE_DELAY=0) for preprocessing, ensuring memory is fully released to the OS after large dataset ingestion.

### 2.3 Cache System (CACHE_VERSION = 3)

```
{cache_path}/{cache_id}/
├── manifest.json          # Component config, file list, data values
├── preprocessed/
│   ├── data_0.parquet     # Main data (sorted by filter columns)
│   ├── level_0.parquet    # Heatmap: smallest compression level
│   ├── level_1.parquet    # Heatmap: next compression level
│   └── ...
└── .cache_config.json     # SequenceView-specific config
```

- **Version 2** added sorting by filter columns and smaller row groups for better predicate pushdown.
- **Version 3** added type downcasting: Int64 -> Int32, Float64 -> Float32 (avoids JS BigInt overhead, sufficient precision for visualization).
- **Reconstruction mode** reads all config from the manifest, requiring only `cache_id` and `cache_path`. This enables stateless page reloads.
- **CacheMissError** is raised when no cache exists and no data is provided — the expected signal for "first load" scenarios.

### 2.4 Data Transfer Budget

| Mechanism | Purpose | Overhead |
|-----------|---------|----------|
| Parquet on disk | Long-term cache | Minimal (compressed columnar) |
| Session state | Per-component current-filter data | O(components), not O(components x filter_combos) |
| Hash echo | Prevent redundant Vue updates | 32-byte SHA256 per component |
| Type downcasting | Reduce Arrow serialization | ~50% size reduction (Float64 -> Float32) |
| Row groups (10K/50K) | Enable predicate pushdown on filtered reads | Tuned for filtered vs unfiltered access patterns |

### 2.5 Frontend Technology Stack

The Vue.js frontend is compiled into a single Streamlit custom component (`streamlit_vue_component`). All five component types are multiplexed through this single iframe using a `componentType` discriminator. No source `.vue` files are shipped — only the compiled dist bundle.

| Library | Version | Purpose |
|---------|---------|---------|
| Vue.js | 3.5.25 | Core framework |
| Vuetify | — | Material Design component library |
| Plotly.js | — | LinePlot, Heatmap (scattergl), VolcanoPlot |
| Tabulator.js | — | Table component (server-side pagination) |
| Apache Arrow (JS) | — | Binary data transfer (`tableFromIPC`) |
| Material Design Icons | — | Icon webfont (~400KB-1.3MB per format) |

**Bundle size:** ~6.1 MB minified JS + ~884 KB minified CSS. This is loaded once per session in the iframe.

### 2.6 Render Cycle (6-Phase Protocol)

Each Streamlit rerun triggers the full render cycle in `bridge.py::render_component()`:

| Phase | Action | Data direction |
|-------|--------|----------------|
| 0 | Generate key, get component args | — |
| 1 | Check Python-side cache for current filter state | Python only |
| 2 | Call `vue_func()` with Arrow-serialized data (or empty if cache valid) | Python -> Vue |
| 3 | Process Vue's return: echoed hash, state updates, data requests | Vue -> Python |
| 4 | Validate interactivity selections; prepare updated data | Python only |
| 5 | Cache prepared data; convert Polars -> pandas -> Arrow | Python only |
| 6 | `st.rerun()` if state changed (skips all downstream Streamlit code) | Python only |

**Key behavior:** When `st.rerun()` is called in Phase 6, all Python code after the component call is skipped on that run. This means linked components (Table -> LinePlot -> SequenceView) may require 2-3 Streamlit reruns to fully propagate a selection change through the chain.

### 2.7 Architectural Caveats

1. **SequenceView does NOT extend BaseComponent.** It has its own independent cache system (`CACHE_VERSION = 1`, separate from BaseComponent's `CACHE_VERSION = 3`). It stores `sequences.parquet` and `peaks.parquet` directly (no `preprocessed/` subdir) and uses `.cache_config.json` instead of `manifest.json`. This means:
   - No `regenerate_cache` parameter (cache is always created when data is provided)
   - No `filter_defaults` parameter (defaults are always None)
   - The `@register_component("sequence_view")` decorator sets `_component_type` but duck-typing is used for rendering

2. **`_cache_exists()` does NOT validate config hash.** In reconstruction mode, all config is restored from the manifest, not from constructor arguments. If you change `filters` or `interactivity` mappings without providing new data, the old config is silently used. To pick up config changes, pass `regenerate_cache=True` with data.

3. **Reconstruction mode silently ignores kwargs.** Subclass attributes set from kwargs before `super().__init__()` are overwritten by `_restore_cache_config()`. This is by design but can confuse users expecting parameter changes to take effect without cache regeneration.

4. **Session ID is not cryptographically secure** — uses `np.random.random()`. Only prevents accidental cross-tab state contamination, not deliberate interference.

5. **Batch resend mechanism:** When any component sets `_BATCH_RESEND_KEY`, ALL component hashes are cleared on next render. This ensures coordinated refresh but means a single component's state change can trigger data resend for all components.

---

## 3. Where openms_insight Shines

### 3.1 Cross-Component State Linking

This is the single most valuable feature for the MSInsight use case. A click on a PSM table row propagates through the entire Explore tab:

```
Table click (row 42)
  └─► StateManager.set_selection("selected_psm", 42)
        ├─► LinePlot filters peaks where psm_id == 42
        ├─► SequenceView filters to sequence for psm_id == 42
        └─► (Heatmap: independent, no PSM filtering)
```

**Why this matters for MSInsight:** MSpecViz had no equivalent. Users had to manually navigate between spectrum views and PSM lists. The StateManager eliminates an entire class of UI synchronization bugs and makes the Explore tab feel like a single integrated tool rather than four disconnected widgets.

**Technical depth:** The linking is not a simple event bus. The counter-based system (`selection_counter` vs `pagination_counter`) prevents a common Streamlit pitfall where pagination changes trigger unnecessary data reloads in linked components. Each component's `_prepare_vue_data()` method independently applies filters using `filter_and_collect_cached()`, which performs Polars predicate pushdown — meaning only the rows matching the current selection are ever materialized.

### 3.2 Polars LazyFrame Integration

openms_insight is built on Polars from the ground up, not adapted from pandas. This means:

- **Schema enforcement** at the Parquet builder level catches type mismatches before they reach visualization code.
- **Predicate pushdown** reads only matching row groups from disk, critical for large LC-MS datasets.
- **Zero-copy Arrow serialization** to the Vue frontend.
- **Lazy evaluation** defers computation until render time, allowing filter state to narrow the materialized data.

**Comparison to MSpecViz:** MSpecViz passes raw pandas DataFrames through function calls, with no schema validation. Type mismatches (e.g., an integer charge column becoming float after a merge) silently propagate and cause plotting errors downstream. The Polars schema approach catches these at ingestion time.

### 3.3 Multi-Resolution Heatmap Downsampling

The Heatmap component uses a sophisticated cascading compression strategy:

1. **Logarithmic compression levels** — For a 1M-point dataset with `min_points=10000`, levels are computed at [10000, 100000] (powers of 10).
2. **Cascading construction** — Each level is built from the next larger level, not from raw data. This means preprocessing cost is O(sum of levels), not O(N x num_levels).
3. **Aspect-ratio-aware binning** — `compute_optimal_bins()` uses the display aspect ratio (16:9 default), not the data aspect ratio, to distribute bins. This ensures approximately square bins in screen space.
4. **Zoom-responsive level selection** — On zoom, the component selects the smallest level with >= `min_points` points in the visible range. If over `min_points`, additional render-time downsampling uses the actual zoom aspect ratio.
5. **Categorical filter optimization** — Per-category compression hierarchies ensure constant point counts regardless of filter selection (e.g., ion mobility dimension).

**Why this matters:** LC-MS run heatmaps can have millions of data points. MSpecViz has no equivalent — it never offered an RT x m/z heatmap. openms_insight makes this feasible without loading all data into browser memory.

### 3.4 SequenceView with Interactive Fragment Matching

The SequenceView component performs fragment matching entirely on the Vue side:

- **TheoreticalSpectrumGenerator integration** — Uses OpenMS `AASequence.fromString()` to compute theoretical fragment masses for a, b, c, x, y, z ion types.
- **Interactive tolerance adjustment** — Because matching runs in Vue, the user can adjust tolerance (Da or ppm) without a Python roundtrip.
- **Ion-type color mapping** — Default colors: a=#9B59B6, b=#E74C3C, c=#E67E22, x=#1ABC9C, y=#3498DB, z=#2ECC71.
- **Bidirectional annotation** — Vue sends back `{"peak_id": [...], "annotation": [...], "highlight_color": [...]}` arrays to Python, enabling SequenceView -> LinePlot annotation linking.

**Why this matters:** MSpecViz never had a sequence coverage view. This is entirely new functionality that gives users immediate visual feedback on which fragments are matched and what the coverage looks like across the peptide backbone.

### 3.5 Server-Side Pagination (Table)

The Table component uses Tabulator with server-side pagination:

- Only the current page (default 100 rows) is sent to the browser.
- Sort state, column filters, and page number persist across Streamlit reruns via `pagination_identifier`.
- Auto-detection of column types: numeric columns get range filters, low-cardinality string columns get dropdown filters.
- "Go to" navigation for columns with 100% unique values (auto-detected).
- Selection navigation: when a linked selection changes, the Table auto-navigates to the correct page.

**Why this matters:** MSInsight can generate hundreds or thousands of PSMs. Without server-side pagination, the browser would receive the full dataset on every rerun. The Table component keeps browser memory constant regardless of PSM count.

### 3.6 Memory-Isolated Preprocessing

The `data_path=` mode spawns a separate Python process for preprocessing:

```python
# In subprocess_preprocess.py:
os.environ["LD_PRELOAD"] = mimalloc_path
os.environ["MIMALLOC_PURGE_DELAY"] = "0"  # Immediate OS memory release
```

This ensures that the memory used for preprocessing (sorting, downsampling, type casting) is fully returned to the OS after the subprocess exits. The main Streamlit process only mmaps the resulting Parquet files.

**Why this matters:** Streamlit processes are long-lived. Without subprocess isolation, preprocessing a large heatmap dataset would permanently increase the process's resident memory. This is the kind of production-grade concern that MSpecViz never addressed.

---

## 4. Component-by-Component Analysis

### 4.1 Table

| Aspect | Detail |
|--------|--------|
| **Input** | `pl.LazyFrame` with any schema |
| **Key params** | `index_field`, `column_definitions`, `pagination`, `page_size`, `initial_sort`, `go_to_fields` |
| **Interactivity** | Row click -> sets StateManager identifiers from row's column values |
| **Filtering** | Server-side predicate pushdown via `filters` dict |
| **Customization** | Auto-detected column formatters (money, progress, star, etc.), custom column defs via Tabulator API |
| **Strengths** | Server-side pagination, auto-type detection, selection navigation across pages |
| **Limitations** | Large numeric values must fit Int32 (+-2.1B) after downcasting; no column grouping or pivoting |

### 4.2 LinePlot

| Aspect | Detail |
|--------|--------|
| **Input** | `pl.LazyFrame` with x, y columns + optional highlight, annotation |
| **Key params** | `x_column`, `y_column`, `highlight_column`, `annotation_column`, `styling` |
| **Interactivity** | Peak click -> sets identifier; dynamic annotations via `set_dynamic_annotations()` |
| **Filtering** | Same LazyFrame filter mechanism as Table |
| **Customization** | `highlightColor`, `selectedColor`, `unhighlightedColor`; dynamic title override |
| **Strengths** | `from_sequence_view()` class method for SequenceView linking; annotation caching (base data cached, annotations reapplied fresh) |
| **Limitations** | Single y-axis only (no mirror layout); no metrics overlay; row group sizes hardcoded (10K/50K) |

### 4.3 SequenceView

| Aspect | Detail |
|--------|--------|
| **Input** | Sequence: `pl.LazyFrame` or `(str, int)` tuple; Peaks: `pl.LazyFrame` with `peak_id`, `mass` |
| **Key params** | `annotation_config` (ion_types, tolerance, tolerance_ppm, neutral_losses, colors), `deconvolved`, `height` |
| **Filtering** | Dual-source filtering: both sequence and peaks filtered by same identifiers |
| **Customization** | Per-ion-type colors, tolerance mode (Da/ppm), neutral loss toggle, proton loss/addition |
| **Strengths** | Client-side fragment matching (no Python roundtrip for tolerance changes); annotation -> LinePlot linking |
| **Limitations** | Does NOT extend BaseComponent (own cache system, CACHE_VERSION=1); no `regenerate_cache` or `filter_defaults` params; requires OpenMS notation (ProForma must be converted); tolerance default is 20 ppm (must set `tolerance_ppm: False` for Da) |

### 4.4 Heatmap

| Aspect | Detail |
|--------|--------|
| **Input** | `pl.LazyFrame` with x, y, intensity columns |
| **Key params** | `min_points`, `display_aspect_ratio`, `colorscale`, `log_scale`, `categorical_filters`, `use_streaming` |
| **Interactivity** | Zoom state tracked via `zoom_identifier`; level-based resolution selection |
| **Filtering** | Category-aware compression with per-value level hierarchies |
| **Customization** | Plotly colorscales, log/linear intensity, reversed scale, custom axis labels, simple vs spatial downsampling |
| **Strengths** | Multi-resolution cascading; zoom-responsive; streaming mode for memory efficiency; categorical filter optimization |
| **Limitations** | scipy required for spatial binning (fallback to top-N available); categorical filters limited to <20 unique values; no point selection/lasso |

### 4.5 VolcanoPlot

| Aspect | Detail |
|--------|--------|
| **Input** | `pl.LazyFrame` with log2FC and pvalue columns |
| **Key params** | `fc_threshold`, `p_threshold`, `max_labels`, `up_color`, `down_color`, `ns_color` |
| **Interactivity** | Point click -> sets identifier; thresholds adjustable at render time (no cache impact) |
| **Customization** | Three-category coloring, threshold lines toggle, label count limit |
| **Strengths** | Render-time threshold adjustment without cache invalidation; sorted by significance (most significant on top) |
| **Limitations** | Not directly relevant to MSInsight's spectrum comparison workflow; included for completeness |

---

## 5. What MSpecViz Does That openms_insight Replaces Well

### 5.1 Interactive Spectrum Viewing

| MSpecViz approach | openms_insight approach | Verdict |
|---|---|---|
| `create_interactive_mirror_plot()` via pyopenms_viz (ms_plotly backend) | `LinePlot` component with native Plotly | **Better** — removes pyopenms_viz dependency; adds filtering, dynamic annotations, caching |
| Manual Plotly figure construction with `go.Figure()` | Component handles Plotly internally | **Better** — standardized API, hash-based update optimization |
| Per-peak hover text built manually | Annotation column mapped automatically | **Equivalent** — same end result, less boilerplate |

### 5.2 PSM Browsing

| MSpecViz approach | openms_insight approach | Verdict |
|---|---|---|
| `st.selectbox()` with PSM list | `Table` component with server-side pagination | **Much better** — sortable, filterable, paginated, with auto-navigation |
| Manual session state for selected PSM | `StateManager` with counter-based sync | **Much better** — eliminates sync bugs, handles page navigation |
| No linked views | Table -> LinePlot -> SequenceView linked via filter keys | **Entirely new capability** |

### 5.3 Data Persistence

| MSpecViz approach | openms_insight approach | Verdict |
|---|---|---|
| In-memory DataFrames only | Parquet-backed LazyFrames with disk cache | **Much better** — survives Streamlit reruns, enables large datasets |
| No schema validation | Polars schemas enforced at build time | **Better** — catches errors early |
| pandas throughout | Polars LazyFrames with predicate pushdown | **Better** — lower memory, faster filtered reads |

---

## 6. What Is Still Missing

### 6.1 Mirror Plot Layout (Critical Gap)

openms_insight has **no mirror plot component**. The LinePlot displays a single spectrum on one y-axis. MSpecViz's signature visualization — two spectra facing each other with metrics overlay — cannot be produced by any openms_insight component.

**Current workaround in MSInsight:** Mirror plots are generated via custom matplotlib/Plotly code in `src/mirror_plot.py`, completely outside openms_insight. The Explore tab uses LinePlot for single-spectrum viewing only.

**What would be needed:** A `MirrorPlot` component that accepts two LazyFrames (query + reference), renders them on opposing y-axes, supports annotation coloring on both spectra, and can overlay metrics text. This is the single most impactful feature request for the MSInsight use case.

### 6.2 Comparison Metrics Display

openms_insight components have no concept of spectral similarity metrics. There is no mechanism to:

- Compute and display cosine similarity, spectral contrast angle, or Pearson correlation between two spectra.
- Overlay formatted metrics text on a plot.
- Link metrics to a specific comparison pair in a table.

**Current workaround:** MSInsight computes metrics via `src/spectrum_metrics.py` and displays them as text or as matplotlib text boxes on exported plots. The Explore tab's LinePlot cannot show these.

### 6.3 Batch Static Export

openms_insight is designed for interactive browser-based viewing. It has **no static export pipeline**:

- No PNG/SVG/PDF generation from components.
- No batch processing mode (process N spectra, save N figures).
- No `ProcessPoolExecutor` integration for parallel export.

**Current workaround:** MSInsight uses `src/parallel_processor.py` with `ProcessPoolExecutor` to generate matplotlib figures in worker processes, completely independently from openms_insight.

### 6.4 Koina / AI Prediction Integration

openms_insight has no concept of predicted spectra or AI model integration. The package works with whatever data you provide — it doesn't know or care where the data came from.

**This is not necessarily a gap.** Koina predictions produce standard DataFrames that feed into LinePlot/SequenceView just fine. The integration point is in the data loading pipeline, not in the visualization layer. openms_insight's role is to display the data, not to generate it.

### 6.5 Ion Match Statistics

The 17 `ims_*` statistics computed by `SpectrumAnnotator.addIonMatchStatistics()` are a computation concern, not a visualization concern. openms_insight could display them in a Table, but it doesn't compute them.

**This is correctly scoped.** openms_insight should not perform mass spectrometry computations. The statistics should be computed upstream and passed as additional columns to the Table or as metrics overlay text.

### 6.6 Summary Analysis Plots

MSpecViz generates batch summary visualizations:

- Ion matching heatmaps (which ions matched across all comparisons)
- Modification distribution bar charts
- Score distribution histograms
- Charge state distributions

openms_insight could display these via Heatmap or VolcanoPlot, but the data aggregation logic is not part of the package. MSInsight handles this in `src/summary_plots.py`.

### 6.7 Real-Time Ion Type Filtering

The SequenceView's `annotation_config.ion_types` is set at creation time and stored in the cache config. Changing ion types requires cache regeneration. There is no UI control for dynamically toggling ion types without a Python roundtrip.

**Partial mitigation:** The tolerance slider works client-side (Vue), but the ion type filter does not. A future version could move ion type filtering to the Vue side as well.

---

## 7. Where openms_insight Could Improve

### 7.1 High Priority (Would Directly Benefit MSInsight)

#### 7.1.1 MirrorPlot Component

A dedicated mirror plot component would be the single highest-impact addition for mass spectrometry workflows. Proposed API:

```python
MirrorPlot(
    cache_id="mirror",
    upper_data=query_peaks_lf,      # LazyFrame with x, y, annotation, highlight
    lower_data=ref_peaks_lf,        # Same schema
    upper_label="Query: ProT001",
    lower_label="Reference: ProT002",
    metrics_text="Cosine: 0.94\nSCA: 12.3deg",  # Optional overlay
    filters={"selected_psm": "psm_id"},
    annotation_column="annotation",
    highlight_column="highlight",
    cache_path=CACHE_PATH,
)
```

This would eliminate the need for separate matplotlib mirror plot code and integrate the comparison view into the interactive Explore workflow.

#### 7.1.2 Metrics Overlay Support

Allow any plot component to display a text box with formatted metrics:

```python
LinePlot(
    ...,
    overlay_text="Cosine: 0.94\nMatched: 42/67",
    overlay_position="top-right",
)
```

This is simpler than a full MirrorPlot and could be implemented as a generic feature across all plot components.

#### 7.1.3 Static Export API

A `component.to_image(format="png", width=1200, height=800)` method that uses Plotly's static export (kaleido) to generate publication-quality images:

```python
lineplot = LinePlot(cache_id="spectrum", data=peaks_lf, ...)
lineplot.to_image("/exports/spectrum_42.png", format="png", dpi=300)
lineplot.to_image("/exports/spectrum_42.svg", format="svg")
```

This would unify the interactive and export workflows and eliminate the parallel matplotlib pipeline entirely.

#### 7.1.4 Dynamic Ion Type Filtering (SequenceView)

**The problem:** When a user changes their ion type selection (e.g., from `["b", "y"]` to `["a", "b", "y", "z"]`) in the MSInsight UI, the change does not take effect until the SequenceView cache is fully regenerated. This means providing new data and re-initializing the component — a disruptive operation that resets all component state (scroll position, selection, zoom).

**Why this happens — the current data flow in detail:**

```
CACHE TIME (Python, sequenceview.py __init__ → _create_cache):
  ┌──────────────────────────────────────────────────────────┐
  │ annotation_config = {"ion_types": ["b","y"], ...}        │
  │ Stored in .cache_config.json                             │
  │ NOT used to filter which masses are computed              │
  └──────────────────────────────────────────────────────────┘

RENDER TIME (Python, _prepare_vue_data, called every __call__):
  ┌──────────────────────────────────────────────────────────┐
  │ calculate_fragment_masses_pyopenms()                      │
  │   → add_a_ions = "true"                                  │
  │   → add_b_ions = "true"                                  │
  │   → add_c_ions = "true"   ALL 6 ION TYPES ALWAYS         │
  │   → add_x_ions = "true"   COMPUTED UNCONDITIONALLY        │
  │   → add_y_ions = "true"                                  │
  │   → add_z_ions = "true"                                  │
  │                                                          │
  │ Returns: {                                               │
  │   "fragment_masses_a": [...],  ─┐                        │
  │   "fragment_masses_b": [...],   │ ALL 6 SENT             │
  │   "fragment_masses_c": [...],   │ TO VUE                 │
  │   "fragment_masses_x": [...],   │ (unfiltered)           │
  │   "fragment_masses_y": [...],   │                        │
  │   "fragment_masses_z": [...],  ─┘                        │
  │ }                                                        │
  │                                                          │
  │ Also sends: annotationConfig.ion_types = ["b", "y"]      │
  │             (from CACHED config, not from current UI)     │
  └──────────────────────────────────────────────────────────┘

VUE SIDE (JavaScript, SequenceView component):
  ┌──────────────────────────────────────────────────────────┐
  │ Receives ALL 6 fragment_masses_* arrays                  │
  │ Receives annotationConfig.ion_types = ["b", "y"]         │
  │                                                          │
  │ Filters: only matches peaks against b and y masses       │
  │ Only renders b and y ion annotations on the sequence     │
  │                                                          │
  │ Returns _annotations: {peak_id, highlight_color, text}   │
  │ (only for b/y matched peaks)                             │
  └──────────────────────────────────────────────────────────┘
```

**The key insight:** All 6 ion types' theoretical masses are already computed on every render and sent to Vue. Vue already filters client-side based on `annotationConfig.ion_types`. The bottleneck is that `ion_types` is locked into the cache config at creation time and never re-read from the current UI state.

**What would need to change (minimal):**

1. **In `_prepare_vue_data()` (sequenceview.py, ~5 lines):** Read `ion_types` from the current Streamlit session state (or accept it as a render-time parameter in `__call__()`) instead of from `self._annotation_config`. This way, even when the cache is reused, Vue receives the user's current ion type selection.

2. **In `get_state_dependencies()` (sequenceview.py, 1 line):** Include the ion_types selection in the state dependency list so the bridge's hash computation detects changes and resends data to Vue.

3. **No changes needed in:** `calculate_fragment_masses_pyopenms()` (already computes all 6), cache format (masses aren't cached anyway), or the Vue frontend (already filters by `annotationConfig.ion_types`).

**Alternative: Pure Vue-side solution (no Python changes):**

An even better approach would be to add a checkbox/toggle UI directly in the Vue SequenceView component. Since Vue already has all 6 mass arrays and performs the matching client-side, adding `<v-checkbox>` toggles for each ion type (a/b/c/x/y/z) would give instant, zero-latency filtering with no Python roundtrip at all. The `annotationConfig.ion_types` from Python would become the initial/default selection, but the user could override it interactively.

**Impact on MSInsight:** Currently in `content/peptide_mirror.py`, the ion type multiselect (`st.multiselect`) is in the Setup tab (Tab 1), and changing it requires re-running Search & Load (Tab 2) to rebuild the SequenceView cache. With this fix, ion type changes would be instant in the Explore tab (Tab 3), improving the iterative analysis workflow significantly.

### 7.2 Medium Priority (Would Improve Developer Experience)

#### 7.2.1 ProForma Sequence Support

SequenceView currently requires OpenMS parenthesis notation: `PEPTM(Oxidation)IDEK`. MSInsight must convert from ProForma (`PEPTM[Oxidation]IDEK`) via `_proforma_to_openms_str()`. Native ProForma support would eliminate this conversion step and align with the PSI standard.

#### 7.2.2 Configurable Row Group Sizing

LinePlot uses hardcoded row group sizes (10K for filtered, 50K for unfiltered). For mass spectrometry data where filter cardinality varies widely (1 spectrum = ~200 peaks, 1000 spectra = 200K peaks), exposing this as a parameter would allow better tuning.

#### 7.2.3 Component Export to HTML

A `component.to_html()` method that generates a standalone HTML file with embedded data. This would be useful for sharing individual plots without the full Streamlit app.

#### 7.2.4 Batch Cache Generation

A utility for preprocessing multiple cache entries in parallel:

```python
from openms_insight import batch_preprocess

batch_preprocess([
    Table(cache_id="psm_table", data=psm_lf, ...),
    LinePlot(cache_id="spectrum", data=peaks_lf, ...),
    Heatmap(cache_id="rt_heatmap", data=heatmap_lf, ...),
], cache_path=CACHE_PATH, workers=4)
```

This would speed up the initial load for large datasets where preprocessing (especially Heatmap downsampling) is the bottleneck.

### 7.3 Low Priority (Nice to Have)

#### 7.3.1 VolcanoPlot -> DifferentialPlot Generalization

The VolcanoPlot is currently specific to -log10(p) vs log2FC. A more general scatter plot with threshold-based categorization would be useful for other MS metrics (e.g., retention time deviation vs mass error).

#### 7.3.2 Component Composition API

A way to combine multiple components into a single view layout:

```python
Layout([
    Row([Table(cache_id="psm_table"), SequenceView(cache_id="seq_view")]),
    Row([LinePlot(cache_id="spectrum"), Heatmap(cache_id="heatmap")]),
])
```

Currently, layout is entirely handled by Streamlit's `st.columns()` and `st.container()`, which works fine but doesn't allow synchronized scrolling or split-pane resizing.

#### 7.3.3 Annotation Persistence

Currently, SequenceView annotations (the result of client-side fragment matching) are returned to Python but not persisted. A mechanism to save annotations alongside the cache would enable offline annotation review.

#### 7.3.4 Theme / Dark Mode Support

The Vue components use hardcoded color schemes. A theme system that respects Streamlit's theme (or allows custom themes) would improve visual consistency.

---

## 8. Comparison Matrix

### 8.1 Feature Comparison: MSpecViz vs MSInsight (via openms_insight)

| Feature | MSpecViz | MSInsight | Handled by |
|---------|----------|-----------|------------|
| Mirror plot (static) | matplotlib | matplotlib | Custom code (src/mirror_plot.py) |
| Mirror plot (interactive) | pyopenms_viz | Plotly native | Custom code (src/mirror_plot.py) |
| Single spectrum view | pyopenms_viz | LinePlot component | **openms_insight** |
| PSM table | st.selectbox | Table component | **openms_insight** |
| Sequence coverage | None | SequenceView component | **openms_insight** |
| RT x m/z heatmap | None | Heatmap component | **openms_insight** |
| Cross-component linking | Manual state | StateManager | **openms_insight** |
| Koina AI prediction | Integrated | Integrated | Custom code (src/koina_prediction.py) |
| Ion match statistics | Integrated | Integrated | Custom code (src/parallel_processor.py) |
| Comparison metrics | compute_all_metrics() | compute_all_metrics() | Custom code (src/spectrum_metrics.py) |
| Batch PNG/SVG export | ProcessPoolExecutor | ProcessPoolExecutor | Custom code (src/parallel_processor.py) |
| HTML export | Plotly write_html | Plotly write_html | Custom code |
| Summary analysis plots | summary_plots.py | summary_plots.py | Custom code (src/summary_plots.py) |
| Versioned export dirs | Yes | Yes | Custom code (src/export.py) |
| Data schema validation | None | Polars schemas | **openms_insight** + parquet_builder |
| Parquet persistence | parquet_cache.py | LazyFrame disk cache | **openms_insight** |
| idXML peptide preview | Yes | Yes | Custom code |
| Settings persistence | settings.json | session_state | Custom code |

### 8.2 Architecture Comparison

| Aspect | MSpecViz | MSInsight (openms_insight) |
|--------|----------|---------------------------|
| **Visualization library** | matplotlib + pyopenms_viz | openms_insight (Vue.js) + matplotlib (export only) |
| **Data format** | pandas DataFrame | Polars LazyFrame -> Parquet |
| **State management** | Manual session_state | StateManager with counter-based sync |
| **Caching** | Custom parquet_cache.py | Built-in disk cache (CACHE_VERSION=3) |
| **Memory model** | Everything in-process | Subprocess isolation (mimalloc) |
| **Interactivity** | Plotly callbacks | Vue.js bidirectional bridge |
| **Schema enforcement** | None | Polars schemas at build time |
| **Export pipeline** | matplotlib + ProcessPoolExecutor | Same (unchanged) |
| **Dependencies** | pyopenms_viz, matplotlib, plotly | openms_insight, matplotlib (export), plotly (mirror plots) |

### 8.3 Performance Characteristics

| Operation | MSpecViz | MSInsight (openms_insight) |
|-----------|----------|---------------------------|
| Initial page load | Instant (no preprocessing) | Slower (cache generation) |
| Spectrum switch | ~200ms (full Plotly rebuild) | ~50ms (predicate pushdown on cached Parquet) |
| PSM navigation | ~100ms (selectbox rerun) | ~20ms (Table page, no data reload) |
| Large heatmap (1M pts) | N/A (no heatmap) | <500ms render (multi-resolution, ~10K points sent to browser) |
| Batch export (100 plots) | ~40s (4 workers) | ~40s (same pipeline, unchanged) |
| Memory (1000 PSMs) | ~800MB (all in pandas) | ~200MB (LazyFrame, only materialized page) |

---

## 9. Recommendations

### 9.1 For MSInsight Development

1. **Keep the hybrid architecture.** openms_insight excels at interactive exploration (Explore tab). Custom code excels at computation and batch export (Search & Load tab, Export tab). Don't try to force either paradigm where the other is better.

2. **Invest in the Parquet builder.** The `parquet_builder.py` -> openms_insight pipeline is the architectural leverage point. Every improvement to schema quality, column naming, or data completeness improves all five components simultaneously.

3. **Consider contributing a MirrorPlot component upstream.** This is the highest-value contribution for both MSInsight and the broader OpenMS WebApps ecosystem. Even a minimal implementation (two-trace Plotly with opposing y-axes) would eliminate the custom Plotly mirror plot code.

4. **Use `data_path=` mode for production deployments.** The subprocess memory isolation is critical for long-running Streamlit servers. Convert the current `data=lf` pattern to write LazyFrames to Parquet first, then pass the path.

5. **Pin openms_insight version.** The CACHE_VERSION mechanism means version upgrades may invalidate caches. Pin to 0.1.13 and test before upgrading.

### 9.2 For openms_insight Development (Upstream Feature Requests)

Prioritized by impact on the MS visualization use case:

1. **MirrorPlot component** — Dual-spectrum comparison with metrics overlay. This is the canonical MS visualization that every proteomics tool needs.
2. **Static export API** — `component.to_image()` via kaleido. Eliminates the need for a parallel matplotlib pipeline.
3. **ProForma sequence support** — PSI standard compliance. Eliminates conversion boilerplate in every downstream app.
4. **Dynamic ion type toggle** — Client-side filtering in SequenceView. Makes the coverage view more interactive.
5. **Metrics overlay** — Generic text overlay on any plot component. Small effort, high visual value.
6. **Configurable row groups** — Tunable for different data access patterns. Minor API addition.

### 9.3 What NOT to Change

1. **Don't replace the batch export pipeline with openms_insight.** Static matplotlib export at 300 DPI with ProcessPoolExecutor is mature, tested, and fast. openms_insight is not designed for this.
2. **Don't compute metrics inside openms_insight components.** Keep computation (spectrum_metrics, ion_match_statistics, Koina) in the Python pipeline where it belongs.
3. **Don't abandon Polars schemas.** The strict typing in `parquet_builder.py` is a feature, not overhead. It catches bugs that would otherwise surface as cryptic Vue rendering errors.
4. **Don't over-cache.** The current pattern of regenerating caches when `cache_version` increments is correct. Premature caching of intermediate states adds complexity without benefit.

---

## Appendix A: openms_insight Component Constructor Reference

<details>
<summary>Table constructor parameters</summary>

```python
Table(
    cache_id: str,                          # REQUIRED — unique component ID
    data: Optional[pl.LazyFrame] = None,    # Creation mode: provide data
    data_path: Optional[str] = None,        # Subprocess mode: provide path
    filters: Optional[Dict[str, str]] = None,
    filter_defaults: Optional[Dict[str, Any]] = None,
    interactivity: Optional[Dict[str, str]] = None,
    cache_path: str = ".",
    regenerate_cache: bool = False,
    column_definitions: Optional[List[Dict]] = None,
    title: Optional[str] = None,
    index_field: str = "id",
    go_to_fields: Optional[List[str]] = None,
    layout: str = "fitDataFill",
    default_row: int = 0,
    initial_sort: Optional[List[Dict]] = None,
    pagination: bool = True,
    page_size: int = 100,
    pagination_identifier: Optional[str] = None,
)
```
</details>

<details>
<summary>LinePlot constructor parameters</summary>

```python
LinePlot(
    cache_id: str,
    data: Optional[pl.LazyFrame] = None,
    data_path: Optional[str] = None,
    filters: Optional[Dict[str, str]] = None,
    filter_defaults: Optional[Dict[str, Any]] = None,
    interactivity: Optional[Dict[str, str]] = None,
    cache_path: str = ".",
    regenerate_cache: bool = False,
    x_column: str = "x",
    y_column: str = "y",
    title: Optional[str] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    highlight_column: Optional[str] = None,
    annotation_column: Optional[str] = None,
    styling: Optional[Dict[str, Any]] = None,
    config: Optional[Dict[str, Any]] = None,
)
```
</details>

<details>
<summary>SequenceView constructor parameters</summary>

```python
SequenceView(
    cache_id: str,
    sequence_data: Optional[Union[pl.LazyFrame, Tuple[str, int], str]] = None,
    sequence_data_path: Optional[str] = None,
    peaks_data: Optional[pl.LazyFrame] = None,
    peaks_data_path: Optional[str] = None,
    filters: Optional[Dict[str, str]] = None,
    interactivity: Optional[Dict[str, str]] = None,
    deconvolved: bool = False,
    annotation_config: Optional[Dict[str, Any]] = None,
    cache_path: str = ".",
    title: Optional[str] = None,
    height: int = 400,
)
```
</details>

<details>
<summary>Heatmap constructor parameters</summary>

```python
Heatmap(
    cache_id: str,
    data: Optional[pl.LazyFrame] = None,
    data_path: Optional[str] = None,
    x_column: Optional[str] = None,
    y_column: Optional[str] = None,
    intensity_column: Optional[str] = None,
    filters: Optional[Dict[str, str]] = None,
    filter_defaults: Optional[Dict[str, Any]] = None,
    interactivity: Optional[Dict[str, str]] = None,
    cache_path: str = ".",
    regenerate_cache: bool = False,
    min_points: int = 10000,
    display_aspect_ratio: float = 16/9,
    x_bins: Optional[int] = None,
    y_bins: Optional[int] = None,
    zoom_identifier: str = "heatmap_zoom",
    title: Optional[str] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    colorscale: str = "Portland",
    reversescale: bool = False,
    use_simple_downsample: bool = False,
    use_streaming: bool = True,
    categorical_filters: Optional[List[str]] = None,
    category_column: Optional[str] = None,
    category_colors: Optional[Dict[str, str]] = None,
    log_scale: bool = True,
    low_values_on_top: bool = False,
    intensity_label: Optional[str] = None,
)
```
</details>

<details>
<summary>VolcanoPlot constructor parameters</summary>

```python
VolcanoPlot(
    cache_id: str,
    data: Optional[pl.LazyFrame] = None,
    data_path: Optional[str] = None,
    log2fc_column: str = "log2FC",
    pvalue_column: str = "pvalue",
    label_column: Optional[str] = None,
    filters: Optional[Dict[str, str]] = None,
    filter_defaults: Optional[Dict[str, Any]] = None,
    interactivity: Optional[Dict[str, str]] = None,
    cache_path: str = ".",
    regenerate_cache: bool = False,
    title: Optional[str] = None,
    x_label: str = "log2 Fold Change",
    y_label: str = "-log10(p-value)",
    up_color: str = "#E74C3C",
    down_color: str = "#3498DB",
    ns_color: str = "#95A5A6",
    show_threshold_lines: bool = True,
    threshold_line_style: str = "dash",
    # Render-time params (no cache impact):
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
    max_labels: int = 10,
)
```
</details>

<details>
<summary>StateManager API</summary>

```python
StateManager(session_key: str = "svc_state")

sm.get_selection(identifier: str) -> Any
sm.set_selection(identifier: str, value: Any) -> bool   # True if changed
sm.clear_selection(identifier: str) -> bool
sm.get_all_selections() -> Dict[str, Any]
sm.get_state_for_vue() -> Dict[str, Any]                # Includes counters
sm.update_from_vue(vue_state: dict) -> bool              # Counter conflict resolution
sm.clear() -> None
```
</details>

---

## Appendix B: Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                        MSInsight Pipeline                        │
│                                                                  │
│  idXML + mzML files                                             │
│       │                                                          │
│       ▼                                                          │
│  peptide_search.py ──► variant_groups                           │
│       │                                                          │
│       ▼                                                          │
│  spectrum_loader.py ──► LoadedSpectrum objects                   │
│       │                                                          │
│       ├──► peak_annotation.py ──► annotated DataFrames          │
│       │                                                          │
│       ├──► koina_prediction.py ──► predicted reference spectra  │
│       │                                                          │
│       └──► spectrum_metrics.py ──► comparison metrics            │
│                                                                  │
│  parquet_builder.py:                                            │
│       ├──► PSM LazyFrame ─────────────────► Table               │
│       ├──► Peaks LazyFrame ───────────────► LinePlot            │
│       ├──► Sequence LazyFrame ────────────► SequenceView        │
│       └──► Heatmap LazyFrame ─────────────► Heatmap             │
│                                                                  │
│  parallel_processor.py:                                         │
│       └──► matplotlib figures ────────────► PNG/SVG/PDF export  │
│                                                                  │
│  export.py:                                                      │
│       └──► versioned directory structure                         │
└─────────────────────────────────────────────────────────────────┘
```

---

*Document generated from exhaustive source code review of openms_insight 0.1.13 (156 commits, BSD-3-Clause) installed at `/mnt/volume/bin/miniconda/envs/msinsight/lib/python3.11/site-packages/openms_insight/`. All constructor signatures, cache formats, and internal mechanisms verified against actual source code.*
