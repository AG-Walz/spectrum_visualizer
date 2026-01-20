import json
from pathlib import Path

notebook_path = Path("specviz_full_pipeline.ipynb")
with open(notebook_path, "r") as f:
    nb = json.load(f)

# 1. Update Imports (Cell 37)
import_cell_idx = -1
for i, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "code":
        source_str = "".join(cell["source"])
        if "# Import core modules" in source_str and "from shared_imports import pd" in source_str:
            import_cell_idx = i
            break

if import_cell_idx != -1:
    print(f"✅ Found import cell at index {import_cell_idx}")
    # Replace the cell content with updated imports
    new_imports = [
        "# Import core modules\n",
        "from shared_imports import pd, oms, Peptidoform\n",
        "\n",
        "# CRITICAL: Mock Streamlit cache decorators BEFORE importing any modules that use them\n",
        "# This prevents \"No runtime found\" warnings when running in Jupyter notebooks\n",
        "import sys\n",
        "from unittest.mock import MagicMock\n",
        "\n",
        "# Create mock Streamlit module if it doesn't exist or override it\n",
        "if 'streamlit' not in sys.modules:\n",
        "    sys.modules['streamlit'] = MagicMock()\n",
        "\n",
        "# Create no-op cache decorators that just return the function unchanged\n",
        "def noop_decorator(func=None, **kwargs):\n",
        "    \"\"\"No-op decorator for Streamlit cache functions in notebook context.\"\"\"\n",
        "    if func is None:\n",
        "        # Called with arguments: @cache_data(ttl=3600)\n",
        "        return lambda f: f\n",
        "    else:\n",
        "        # Called without arguments: @cache_data\n",
        "        return func\n",
        "\n",
        "# Override Streamlit cache decorators with no-ops\n",
        "sys.modules['streamlit'].cache_data = noop_decorator\n",
        "sys.modules['streamlit'].cache_resource = noop_decorator\n",
        "\n",
        "# Now import spectrum_analysis - cache decorators will be no-ops\n",
        "import spectrum_analysis\n",
        "\n",
        "# Reload the module to pick up the bug fix for bytes decoding\n",
        "import importlib\n",
        "importlib.reload(spectrum_analysis)\n",
        "\n",
        "# IMPORTANT: Patch the cached loader to avoid pickling pyOpenMS objects in notebooks\n",
        "# The app uses @st.cache_data which pickles return values, but pyOpenMS objects can't be pickled\n",
        "# This replacement allows the notebook to work without Streamlit's caching\n",
        "def _load_idxml_peptides_no_cache(idxml_path):\n",
        "    protein_ids = []\n",
        "    peptide_ids = []\n",
        "    oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)\n",
        "    return protein_ids, peptide_ids\n",
        "\n",
        "spectrum_analysis._load_idxml_peptides = _load_idxml_peptides_no_cache\n",
        "\n",
        "# Import the EXACT same pipeline functions used by the app\n",
        "from spectrum_analysis import (\n",
        "    filter_idxml_by_selected_peptides,  # PHASE 2: Filter idXMLs by peptide list\n",
        "    load_only_needed_spectra,           # PHASE 3: Load only needed spectra (memory-efficient)\n",
        ")\n",
        "\n",
        "# Import the EXACT same plotting functions used by the app\n",
        "from spectrum_engine import (\n",
        "    create_matplotlib_mirror_plot,      # Create mirror plot figure\n",
        "    export_matplotlib_mirror_plot,      # Export plot to file\n",
        "    convert_spectrum_to_df,             # Convert spectrum to DataFrame\n",
        "    generate_theoretical_spectrum,      # Generate theoretical spectrum\n",
        "    align_spectra,                      # Align experimental and theoretical spectra\n",
        "    annotate_peaks,                     # Annotate matched peaks\n",
        "    extract_annotations,                # Extract annotations for plotting\n",
        "    create_peptide_overview_plot,       # Create summary barplot\n",
        "    export_overview_plot,               # Export summary plot\n",
        "    EXP_COLORS,                         # Color scheme for experimental spectra\n",
        "    PRED_COLORS,                        # Color scheme for predicted spectra\n",
        ")\n",
        "\n",
        "from correlation import spectrum_correlation, new_spectrum_correlation\n",
        "\n",
        "print(\"✅ Modules loaded successfully (Streamlit cache decorators disabled for notebook)\")\n"
    ]
    nb["cells"][import_cell_idx]["source"] = new_imports
else:
    print("❌ Could not find import cell")

# 2. Insert Summary Logic (Before Summary)
summary_cell_idx = -1
for i, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "markdown":
        source_str = "".join(cell["source"])
        if "## Summary" in source_str:
            summary_cell_idx = i
            break

if summary_cell_idx != -1:
    print(f"✅ Found summary cell at index {summary_cell_idx}")
    
    cell_a = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "print(\"📊 Calculating One-vs-All statistics...\")\n",
            "\n",
            "# Initialize results container\n",
            "r2_results = []\n",
            "summary_plot_output_dir = specviz_dir / \"plots\" / project_name / \"summary\"\n",
            "summary_plot_output_dir.mkdir(parents=True, exist_ok=True)\n",
            "\n",
            "# Loop through all filtered peptides\n",
            "for peptide_seq in peptides:\n",
            "    # Get all rows for this peptide\n",
            "    peptide_rows = mapping_df[mapping_df[\"peptide_sequence\"].str.contains(peptide_seq, regex=False)]\n",
            "    \n",
            "    if peptide_rows.empty:\n",
            "        continue\n",
            "        \n",
            "    # Separate into reference and samples\n",
            "    ref_rows = peptide_rows[peptide_rows[\"mzml_file\"].isin(reference_mzml_names)]\n",
            "    sample_rows = peptide_rows[~peptide_rows[\"mzml_file\"].isin(reference_mzml_names)]\n",
            "    \n",
            "    if ref_rows.empty or sample_rows.empty:\n",
            "        continue\n",
            "        \n",
            "    # For summary, we compare EACH sample against the BEST reference\n",
            "    for _, sample_row in sample_rows.iterrows():\n",
            "        ref_row = pick_reference_row(sample_row)\n",
            "        if ref_row is None:\n",
            "            continue\n",
            "            \n",
            "        # Load spectra\n",
            "        sample_spec = load_spectrum_from_experiment(experiments, lookups, sample_row[\"mzml_file\"], sample_row[\"spectrum_id\"])\n",
            "        ref_spec = load_spectrum_from_experiment(experiments, lookups, ref_row[\"mzml_file\"], ref_row[\"spectrum_id\"])\n",
            "        \n",
            "        if sample_spec is None or ref_spec is None:\n",
            "            continue\n",
            "            \n",
            "        # Annotate\n",
            "        # Note: We assume spectra loaded from experiments are raw and need annotation\n",
            "        sample_df = annotate_spectrum_with_theoretical(\n",
            "             sample_spec, sample_row[\"peptide_sequence\"], charge=sample_row.get(\"charge\"), \n",
            "             tol=fragment_mass_tolerance, included_ions=ion_types\n",
            "        )\n",
            "        ref_df = annotate_spectrum_with_theoretical(\n",
            "             ref_spec, ref_row[\"peptide_sequence\"], charge=ref_row.get(\"charge\"), \n",
            "             tol=fragment_mass_tolerance, included_ions=ion_types\n",
            "        )\n",
            "        \n",
            "        # Calculate R2\n",
            "        # Method 1: Matching\n",
            "        sample_ions = [(r.mz, r.intensity, r.annotation) for _, r in sample_df.iterrows() if r.annotation != 'unassigned']\n",
            "        ref_ions = [(r.mz, r.intensity, r.annotation) for _, r in ref_df.iterrows() if r.annotation != 'unassigned']\n",
            "        \n",
            "        r2_matching = None\n",
            "        if sample_ions and ref_ions:\n",
            "             r2_matching, _, _, _ = spectrum_correlation(tuple(sample_ions), tuple(ref_ions))\n",
            "             \n",
            "        # Method 2: Binned\n",
            "        r2_binned = None\n",
            "        if not sample_df.empty and not ref_df.empty:\n",
            "             r2_binned, _, _, _ = new_spectrum_correlation(sample_df, ref_df)\n",
            "             \n",
            "        # Store result\n",
            "        r2_results.append({\n",
            "            \"peptidoform\": peptide_seq,\n",
            "            \"sample_run\": sample_row['run'],\n",
            "            \"reference_run\": ref_row['run'],\n",
            "            \"r2_matching\": r2_matching if r2_matching is not None else 0,\n",
            "            \"r2_binned\": r2_binned if r2_binned is not None else 0\n",
            "        })\n",
            "\n",
            "r2_data = pd.DataFrame(r2_results)\n",
            "print(f\"✅ Calculated statistics for {len(r2_data)} comparisons\")\n"
        ]
    }
    
    cell_b = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "if not r2_data.empty:\n",
            "    print(\"📈 Generating overview summary plot...\")\n",
            "    \n",
            "    # Create overview plot\n",
            "    # Note: create_peptide_overview_plot expects 'r2_matching' or 'r2_binned' column\n",
            "    fig, summary_df = create_peptide_overview_plot(\n",
            "        r2_data=r2_data,\n",
            "        r2_metric='r2_matching', # Default to matching\n",
            "        top_n=None\n",
            "    )\n",
            "    \n",
            "    if fig:\n",
            "        # Save plot\n",
            "        plot_path = export_overview_plot(fig, str(summary_plot_output_dir), \"matching\")\n",
            "        display(fig)\n",
            "        print(f\"💾 Summary plot saved to: {plot_path}\")\n",
            "    else:\n",
            "        print(\"⚠️ Failed to create summary plot\")\n",
            "        \n",
            "    # Export Excel\n",
            "    excel_path = analysis_output_dir / \"summary_r2_data.xlsx\"\n",
            "    r2_data.to_excel(excel_path, index=False)\n",
            "    print(f\"💾 Summary data saved to: {excel_path}\")\n",
            "    \n",
            "    # Display top 5 rows\n",
            "    display(r2_data.head())\n",
            "else:\n",
            "    print(\"⚠️ No R2 data to plot\")\n"
        ]
    }

    markdown_header = {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Step 8: Calculate One-vs-All Summary\n",
            "\n",
            "Calculate R-squared values for all samples against their best reference match and generate summary statistics and plots."
        ]
    }

    # Insert before Summary cell
    nb["cells"].insert(summary_cell_idx, markdown_header)
    nb["cells"].insert(summary_cell_idx + 1, cell_a)
    nb["cells"].insert(summary_cell_idx + 2, cell_b)
    print("✅ Inserted summary cells")
else:
    print("❌ Could not find summary cell")

# Save
with open(notebook_path, "w") as f:
    json.dump(nb, f, indent=1)
print("✅ Notebook patched and saved.")
