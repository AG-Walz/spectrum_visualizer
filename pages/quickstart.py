"""Quickstart guide - Getting started for new users.

Step-by-step workflow overview with instructions. Static content page.
Accessed via: App sidebar → Home → Quickstart
"""
import streamlit as st
from pathlib import Path

# Simple setup
st.set_page_config(
    page_title="SpecViz - Quickstart",
    initial_sidebar_state="expanded"
)

# Load markdown documentation (single source of truth)
docs_path = Path(__file__).parent.parent / "docs" / "QUICKSTART.md"

try:
    with open(docs_path, "r", encoding="utf-8") as f:
        quickstart_content = f.read()
    st.markdown(quickstart_content, unsafe_allow_html=True)
except FileNotFoundError:
    st.error(f"❌ Documentation file not found: {docs_path}")
    st.info("Please ensure docs/QUICKSTART.md exists in the project root.")

    st.markdown("---")
    
    st.markdown("""
    ## 🚀 Running the App
    
    ```bash
    conda activate specviz
    streamlit run specviz/app.py
    ```
    
    The app will open in your default browser at `http://localhost:8501`.
    
    ---
    ## 🧭 Navigate the App
    
    The sidebar on the left contains all navigation options:
    
    ### 📄 Home Pages
    Static pages with fixed information:
    - **Quickstart** (this page): Quick reference guide to get started
    - <a href="/documentation" target="_self"><b>Documentation</b></a>: Comprehensive documentation with detailed explanations of every feature
    
    ### 🔬 Workspaces
    Interactive environments for analysis (see detailed descriptions below):
    - <a href="/spectrum_analysis" target="_self"><b>Spectrum Analysis</b></a>: Main workspace for MS/MS visualization and comparison
    - <a href="/dataset_filter" target="_self"><b>Dataset Filter</b></a>: Filter and preprocess <a href="/documentation#file-formats" target="_self">mzML files</a>
    - <a href="/settings" target="_self"><b>Settings</b></a>: Configure application preferences (<a href="/documentation#custom-settings" target="_self">custom settings</a>)
    
    **⚠️ Important**: Workspaces maintain state during your session. If you navigate away from a workspace, your progress will be reset. Complete your analysis in one session or export results before leaving.
    
    ---
    
    ## 🔬 Workspace Details
    
    ### Spectrum Analysis
    The primary workspace for visualizing and analyzing MS/MS spectra:
    
    **What it does:**
    - Loads [idXML files](/documentation#file-formats) containing [peptide identifications](/documentation#idxml-processing)
    - Extracts unique [peptidoforms](/documentation#understanding-the-output) and associated spectra
    - Matches spectra to [mzML files](/documentation#mzml-loading)
    - Generates [AI predictions](/documentation#ai-models) using [Koina API](/documentation#koina-api) ([MS2PIP](/documentation#ms2pip), [Prosit](/documentation#prosit), [AlphaPept](/documentation#alphapept))
    - Creates [mirror plots](/documentation#mirror-plot-anatomy) comparing [experimental vs predicted spectra](/documentation#exp-vs-pred-analysis)
    - Performs [multi-spectrum comparisons](/documentation#exp-vs-exp-analysis) (all-vs-all, [one-vs-all](/documentation#one-vs-all-analysis))
    - Exports plots and data to [Excel](/documentation#excel-workbook-contents)/[HTML](/documentation#html-plot-files)
    
    **When to use:**
    - Validating [peptide identifications](/documentation#idxml-processing)
    - Comparing [fragmentation patterns](/documentation#visualization) across [replicates](/documentation#exp-vs-exp-analysis)
    - Evaluating [AI model performance](/documentation#ai-models)
    - [Quality control](/documentation#data-quality-issues) of MS/MS data
    
    **Workflow:** <a href="/documentation#file-upload-section" target="_self">Upload idXML</a> → <a href="/documentation#peptide-selection" target="_self">Select peptides</a> → <a href="/documentation#mzml-loading" target="_self">Filter spectra</a> → <a href="/documentation#ai-model-selection" target="_self">Choose AI model</a> → <a href="/documentation#automated-analysis" target="_self">Generate comparisons</a>
    
    See detailed workflow below and <a href="/documentation" target="_self"><b>full documentation</b></a>.
    
    ### Dataset Filter
    Preprocessing workspace for mzML files:
    
    **What it does:**
    - Loads raw mzML spectrum files
    - Applies quality filters (intensity thresholds, m/z ranges)
    - Removes low-quality or contaminant spectra
    - Exports filtered datasets
    
    **When to use:**
    - Cleaning up [noisy data](/documentation#data-quality-issues) before identification
    - Removing specific m/z ranges (e.g., calibrants)
    - Reducing [file sizes](/documentation#memory-optimization)
    - Preparing data for downstream tools
    
    **📖 Documentation**: See [**File Formats**](/documentation#file-formats) section for [mzML technical specifications](/documentation#file-formats), [**Memory Optimization**](/documentation#memory-optimization), and [**Performance Benchmarks**](/documentation#performance-benchmarks)
    
    ### Settings
    Application configuration workspace:
    
    **What it does:**
    - Adjusts visualization parameters (plot size, colors)
    - Sets default analysis options (mass tolerance, charge states)
    - Configures data directories
    - Manages performance settings
    
    **When to use:**
    - Customizing the [interface](/documentation#visualization) for your workflow
    - Setting up data paths for [batch processing](/documentation#batch-processing)
    - Optimizing for [large datasets](/documentation#memory-issues)
    
    **📖 Documentation**: See [**Advanced Features**](/documentation#advanced-features) for [batch processing](/documentation#batch-processing) and [**custom settings**](/documentation#custom-settings), [**Technical Details**](/documentation#technical-details) for [architecture](/documentation#architecture) and [data flow](/documentation#data-flow), and [**Performance Benchmarks**](/documentation#performance-benchmarks)
    
    ---
    
    ## 🏁 Quick Workflow
    
    Follow these steps to analyze your first spectrum:
    
    ### 1. Navigate to Spectrum Analysis
    - Click **"Workspaces"** in the sidebar
    - Select **"Spectrum Analysis"**
    - You'll see the file upload interface
    
    **📖 Documentation**: See [**File Upload**](/documentation#file-upload-section) section for details, including [**Upload Methods**](/documentation#upload-methods) (browser vs server) and [**Supported File Formats**](/documentation#supported-file-formats)
    
    ### 2. Upload Identification Files
    - Click the **"🔍 Upload Identification Files"** section
    - Select one or more [idXML files](/documentation#file-formats) (hold Shift for multiple)
    - [Supported formats](/documentation#supported-file-formats): `.idxml`, `.mzid`, `.mzidentml`
    - Click **"📋 Process idXML Files"** ([what this does](/documentation#the-process-idxml-files-button))
    - The app will extract all [peptidoforms](/documentation#understanding-the-output) and display [unique sequences](/documentation#peptide-metrics)
    
    **What happens:** The app parses your [identification files](/documentation#file-formats), extracts [peptide sequences with modifications](/documentation#modification-formats), and determines which [mzML files are needed](/documentation#when-mzml-files-are-needed).
    
    **📖 Documentation**: See [**idXML Processing**](/documentation#idxml-processing) section, including [**The Process idXML Files Button**](/documentation#the-process-idxml-files-button) explanation, [**Understanding the Output**](/documentation#understanding-the-output), and [**Modification Formats**](/documentation#modification-formats)
    
    ### 3. Select Peptides
    - A [dropdown](/documentation#the-peptide-dropdown) appears with all unique [peptidoforms](/documentation#understanding-the-output)
    - Search by sequence or [modification](/documentation#modification-formats) (e.g., "M[Oxidation]")
    - Select one or more peptides to analyze ([what happens](/documentation#what-happens-when-you-select))
    - The interface shows [how many spectra match each peptide](/documentation#peptide-metrics)
    
    **Options:**
    - [**Filter mode**](/documentation#filter-modes): Paste a list of peptides to analyze specific sequences
    - [**Browse mode**](/documentation#filter-modes): Explore all available peptides
    
    **📖 Documentation**: See [**Peptide Selection**](/documentation#peptide-selection) section, including [**The Peptide Dropdown**](/documentation#the-peptide-dropdown) interface, [**Filter Modes**](/documentation#filter-modes) (browse vs filter), [**Peptide Metrics**](/documentation#peptide-metrics), and [**Multiple Peptide Analysis**](/documentation#multiple-peptide-analysis) strategies
    
    ### 4. Upload Spectrum Files (if needed)
    - If [mzML files](/documentation#file-formats) aren't on the server, [upload them](/documentation#upload-methods)
    - Only [files matching your selected peptides](/documentation#when-mzml-files-are-needed) are required
    - The app shows [which files are needed](/documentation#file-mapping-validation)
    - Click **"🔄 Process Files"** ([what this does](/documentation#the-process-files-button)) to load spectra
    
    **What happens:** The app loads only the spectra corresponding to your selected peptides, [optimizing memory usage](/documentation#memory-optimization).
    
    **📖 Documentation**: See [**mzML Loading**](/documentation#mzml-loading) section, including [**When mzML Files Are Needed**](/documentation#when-mzml-files-are-needed), [**The Process Files Button**](/documentation#the-process-files-button), [**Understanding the Mapping Table**](/documentation#understanding-the-mapping-table), and [**Troubleshooting File Loading**](/documentation#troubleshooting-file-loading)
    
    ### 5. Configure Analysis
    - **[Choose AI Model](/documentation#ai-model-selection)**: Select from dropdown (recommended: [`ms2pip_timsTOF2024`](/documentation#ms2pip-models))
    - **[Set Mass Tolerance](/documentation#fragment-mass-tolerance)**: Fragment matching tolerance (default: 0.02 Da)
    - **[Max Fragment Charge](/documentation#maximum-fragment-charge)**: Maximum charge for theoretical fragments (default: 2)
    - **[R² Metric](/documentation#r²-metric-selection)**: Choose `matching_ions`, `binned`, or `both` for correlation
    
    **What each option does:**
    - **[Model](/documentation#ai-model-selection)**: Determines [fragmentation prediction algorithm](/documentation#ai-models) and training data
    - **[Mass tolerance](/documentation#fragment-mass-tolerance)**: Controls how closely experimental and theoretical peaks must match
    - **[Fragment charge](/documentation#maximum-fragment-charge)**: Higher values include multiply-charged fragments
    - **[R² metric](/documentation#r²-metric-selection)**: Method for calculating correlation between experimental and predicted spectra
    
    **📖 Documentation**: See [**Analysis Configuration**](/documentation#analysis-config) section including [**AI Model Selection**](/documentation#ai-model-selection), [**Fragment Mass Tolerance**](/documentation#fragment-mass-tolerance), [**Maximum Fragment Charge**](/documentation#maximum-fragment-charge), [**R² Metric Selection**](/documentation#r²-metric-selection), and [**Debug Mode**](/documentation#debug-mode). Also see [**AI Models**](/documentation#ai-models) section with [**MS2PIP Models**](/documentation#ms2pip-models), [**Prosit Models**](/documentation#prosit-models), and [**AlphaPept Models**](/documentation#alphapept-models) comparisons
    
    ### 6. Run Analysis
    - Click **["🚀 Start Automated Analysis"](/documentation#the-start-automated-analysis-button)**
    - The app generates:
      - **[Exp vs Pred plots](/documentation#exp-vs-pred-analysis)**: Your data compared to [AI predictions](/documentation#ai-models)
      - **[Exp vs Exp plots](/documentation#exp-vs-exp-analysis)**: [Comparisons between your spectra](/documentation#exp-vs-exp-analysis) (if multiple)
      - **[One vs All plots](/documentation#one-vs-all-analysis)**: Each spectrum compared to all others
    - View [interactive mirror plots](/documentation#interactive-features) with [ion annotations](/documentation#peak-annotations)
    - [Export results](/documentation#export-options) to [Excel](/documentation#excel-workbook-contents) or [HTML](/documentation#html-plot-files)
    
    **What happens:** The app queries [Koina API](/documentation#koina-api) for predictions, aligns peaks, calculates [R² scores](/documentation#r²-metric-selection), generates [visualizations](/documentation#visualization), and [compiles results](/documentation#export-generation).
    
    **📖 Documentation**: See [**Automated Analysis**](/documentation#automated-analysis) section including [**The Start Automated Analysis Button**](/documentation#the-start-automated-analysis-button), [**Exp vs Pred Analysis**](/documentation#exp-vs-pred-analysis), [**Exp vs Exp Analysis**](/documentation#exp-vs-exp-analysis), [**One vs All Analysis**](/documentation#one-vs-all-analysis), and [**Export Generation**](/documentation#export-generation). Also see [**Visualization**](/documentation#visualization) with [**Mirror Plot Anatomy**](/documentation#mirror-plot-anatomy), [**Interactive Features**](/documentation#interactive-features), and [**Peak Annotations**](/documentation#peak-annotations). For saving, see [**Export Options**](/documentation#export-options) including [**Excel Workbook Contents**](/documentation#excel-workbook-contents) and [**HTML Plot Files**](/documentation#html-plot-files)

    ---
    
    ## 🆘 Quick Troubleshooting
    
    **["No peptidoforms found"](/documentation#processing-errors)**
    - Check that [idXML files](/documentation#file-formats) contain [peptide identifications](/documentation#idxml-processing)
    - Verify [file format is supported](/documentation#supported-file-formats)
    - **📖 Documentation**: See [**Processing Errors**](/documentation#processing-errors) for common file processing errors
    
    **["No matching spectra"](/documentation#file-mapping-validation)**
    - Ensure [mzML files](/documentation#file-formats) match the experiment
    - Check [file names correspond to identifications](/documentation#troubleshooting-file-loading)
    - **📖 Documentation**: See [**File Mapping Validation**](/documentation#file-mapping-validation) and [**Troubleshooting File Loading**](/documentation#troubleshooting-file-loading) for diagnosis
    
    **["Prediction failed"](/documentation#prediction-errors)**
    - Check internet connection ([Koina API](/documentation#koina-api) requires network)
    - Verify [peptide sequence format](/documentation#modification-formats)
    - Try a [different AI model](/documentation#ai-model-selection)
    - **📖 Documentation**: See [**Prediction Errors**](/documentation#prediction-errors) for troubleshooting and [**AI Models**](/documentation#ai-models) section with [**MS2PIP**](/documentation#ms2pip), [**Prosit**](/documentation#prosit), and [**AlphaPept**](/documentation#alphapept) for model-specific issues
    
    **[Other issues](/documentation#troubleshooting):** Enable [**Debug Mode**](/documentation#debug-mode) in Analysis Settings for detailed logs. See full [**Troubleshooting guide**](/documentation#troubleshooting) including [**File Upload Issues**](/documentation#file-upload-issues), [**Visualization Issues**](/documentation#visualization-issues), [**Memory Issues**](/documentation#memory-issues), and [**Data Quality Issues**](/documentation#data-quality-issues).
    
    ---
    
    ## 💡 Tips
    
    - **[Memory](/documentation#memory-optimization)**: Upload [large files from server](/documentation#upload-methods) instead of browser to avoid [timeouts](/documentation#file-upload-issues)
    - **[Performance](/documentation#performance-benchmarks)**: Analyze [fewer peptides](/documentation#multiple-peptide-analysis) at once for faster results
    - **[Quality](/documentation#data-quality-issues)**: Higher [R² scores](/documentation#r²-metric-selection) indicate better [prediction-experiment agreement](/documentation#exp-vs-pred-analysis)
    - **[Export](/documentation#export-options)**: Save [Excel files](/documentation#excel-workbook-contents) for offline analysis and record-keeping
    
    ---
    
    **Ready to dive deeper?** Explore the <a href="/documentation" target="_self"><b>📖 Complete Documentation</b></a> for detailed guides on every feature.
    
    """, unsafe_allow_html=True)