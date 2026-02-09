import json
from pathlib import Path

notebook_path = Path("/mnt/volume/nacho/specviz/spectrum_visualizer/examples/compare_mzml_plotting.ipynb")

def create_code_cell(source_lines):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line + "\n" for line in source_lines]
    }

with open(notebook_path, "r") as f:
    nb = json.load(f)

# --- 1. Update Configuration Cell ---
config_updated = False
for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = "".join(cell["source"])
        if "INPUT_FILES =" in source and "TARGET_SPECTRUM_ID =" in source:
            if "ANALYSIS_BACKEND =" not in source:
                cell["source"].append("\n# Analysis Backend Selection\n")
                cell["source"].append("# Options: \"both\", \"pyopenms\", \"pyteomics\"\n")
                cell["source"].append("ANALYSIS_BACKEND = \"both\"\n")
                print("✅ Updated Configuration Cell with ANALYSIS_BACKEND.")
            else:
                print("ℹ️ Configuration Cell already contains ANALYSIS_BACKEND.")
            config_updated = True
            break

if not config_updated:
    print("❌ Could not find Configuration Cell.")

# --- 2. Create Plotting Functions Cell Content ---
plotting_functions_source = [
    "def plot_spectra(data, source_name, file_name, spectrum_id):",
    "    \"\"\"Plots MS1 and MS2 for a single backend.\"\"\"",
    "    fig, ax = plt.subplots(1, 2, figsize=(15, 5))",
    "    ax = ax.flatten()",
    "    levels = ['MS1', 'MS2']",
    "    colors = {'pyOpenMS': 'blue', 'pyteomics': 'red'}",
    "    color = colors.get(source_name, 'green')",
    "    ",
    "    for i, level in enumerate(levels):",
    "        level_data = data.get(level)",
    "        if level_data:",
    "            mz, inte, sid = level_data",
    "            ax[i].plot(mz, inte, color=color, alpha=0.8, linewidth=1)",
    "            ax[i].set_title(f\"{level} - {source_name} ({sid})\\nFile: {file_name}\")",
    "            ax[i].set_ylabel(\"Intensity\")",
    "            ax[i].set_xlabel(\"m/z\")",
    "        else:",
    "            ax[i].text(0.5, 0.5, f\"{level} Not Found\", ha='center', va='center')",
    "            ax[i].set_title(f\"{level} - {source_name}: Not Found\")",
    "            ",
    "    plt.tight_layout()",
    "    plt.show()",
    "",
    "def plot_spectra_overlay(data_openms, data_pyteomics, file_name):",
    "    \"\"\"Plots MS1 and MS2 for both backends (2x2 comparison grid).\"\"\"",
    "    fig, ax = plt.subplots(2, 2, figsize=(15, 10))",
    "    ax = ax.flatten()",
    "    ",
    "    plot_defs = [",
    "        {'level': 'MS1', 'source': 'pyOpenMS', 'color': 'blue', 'data': data_openms, 'idx': 0},",
    "        {'level': 'MS1', 'source': 'pyteomics', 'color': 'red', 'data': data_pyteomics, 'idx': 1},",
    "        {'level': 'MS2', 'source': 'pyOpenMS', 'color': 'blue', 'data': data_openms, 'idx': 2},",
    "        {'level': 'MS2', 'source': 'pyteomics', 'color': 'red', 'data': data_pyteomics, 'idx': 3},",
    "    ]",
    "    ",
    "    for p in plot_defs:",
    "        i = p['idx']",
    "        level = p['level']",
    "        source = p['source']",
    "        d = p['data'].get(level)",
    "        ",
    "        if d:",
    "            mz, inte, sid = d",
    "            ax[i].plot(mz, inte, color=p['color'], alpha=0.8, linewidth=1)",
    "            ax[i].set_title(f\"{level} - {source} ({sid})\\nFile: {file_name}\")",
    "            ax[i].set_ylabel(\"Intensity\")",
    "            if i >= 2: ax[i].set_xlabel(\"m/z\")",
    "        else:",
    "            ax[i].text(0.5, 0.5, f\"{level} Not Found\", ha='center', va='center')",
    "            ax[i].set_title(f\"{level} - {source}: Not Found\")",
    "            ",
    "    plt.tight_layout()",
    "    plt.show()"
]

# --- 3. Create Main Loop Cell Content ---
main_loop_source = [
    "import os",
    "import matplotlib.pyplot as plt",
    "import numpy as np # Ensure numpy is available",
    "",
    "# Check if ANALYSIS_BACKEND is defined, default to both if not",
    "if 'ANALYSIS_BACKEND' not in locals():",
    "    ANALYSIS_BACKEND = 'both'",
    "",
    "print(f\"Running Analysis with Backend: {ANALYSIS_BACKEND.upper()}\")",
    "",
    "for mzml_file in INPUT_FILES:",
    "    print(f\"\\n{'='*40}\\nProcessing {Path(mzml_file).name}...\")",
    "    ",
    "    if not os.path.exists(mzml_file):",
    "        print(\"File not found, skipping.\")",
    "        continue",
    "        ",
    "    data_openms = {}",
    "    data_pyteomics = {}",
    "    file_name = Path(mzml_file).name",
    "    ",
    "    # --- LOAD DATA ---",
    "    try:",
    "        if ANALYSIS_BACKEND.lower() in ['both', 'pyopenms']:",
    "            data_openms = get_data_pyopenms(mzml_file, EFFECTIVE_SPECTRUM_ID)",
    "            ",
    "        if ANALYSIS_BACKEND.lower() in ['both', 'pyteomics']:",
    "            data_pyteomics = get_data_pyteomics(mzml_file, EFFECTIVE_SPECTRUM_ID)",
    "            ",
    "    except Exception as e:",
    "        print(f\"⚠️ Error loading data: {e}\")",
    "        continue",
    "        ",
    "    # --- PLOT DATA ---",
    "    if ANALYSIS_BACKEND.lower() == 'pyopenms':",
    "        if data_openms:",
    "            plot_spectra(data_openms, 'pyOpenMS', file_name, EFFECTIVE_SPECTRUM_ID)",
    "        else:",
    "            print(\"No pyOpenMS data loaded.\")",
    "            ",
    "    elif ANALYSIS_BACKEND.lower() == 'pyteomics':",
    "        if data_pyteomics:",
    "            plot_spectra(data_pyteomics, 'pyteomics', file_name, EFFECTIVE_SPECTRUM_ID)",
    "        else:",
    "            print(\"No pyteomics data loaded.\")",
    "            ",
    "    elif ANALYSIS_BACKEND.lower() == 'both':",
    "        if data_openms or data_pyteomics:",
    "            # Comparison stats",
    "            print(\"--- Comparison Results ---\")",
    "            for level in ['MS1', 'MS2']:",
    "                om = data_openms.get(level)",
    "                py = data_pyteomics.get(level)",
    "                if om and py:",
    "                    mz_om, int_om, id_om = om",
    "                    mz_py, int_py, id_py = py",
    "                    match_len = len(mz_om) == len(mz_py)",
    "                    match_val = False",
    "                    if match_len and len(mz_om) > 0:",
    "                        match_val = np.allclose(mz_om, mz_py, atol=1e-4) and np.allclose(int_om, int_py, atol=1e-4)",
    "                    status = \"MATCH\" if match_len and match_val else f\"MISMATCH (Len: {len(mz_om)} vs {len(mz_py)})\"",
    "                    print(f\"{level}: {status} | IDs: {id_om} vs {id_py}\")",
    "                else:",
    "                    print(f\"{level}: Data missing in one or both backends.\")",
    "            ",
    "            plot_spectra_overlay(data_openms, data_pyteomics, file_name)",
    "        else:",
    "            print(\"No data found for either backend.\")"
]

# --- 4. Locate and Replace/Insert Cells ---
main_loop_found = False
new_cells = []

for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = "".join(cell["source"])
        if "for mzml_file in INPUT_FILES:" in source and "plt.subplots" in source:
            # Found the main loop cell.
            # Insert plotting functions cell BEFORE this one
            print("✅ Found Main Loop Cell. Inserting plotting functions and updating loop.")
            new_cells.append(create_code_cell(plotting_functions_source))
            
            # Replace the main loop cell content
            cell["source"] = [line + "\n" for line in main_loop_source]
            cell["outputs"] = [] # Clear outputs
            new_cells.append(cell)
            main_loop_found = True
            continue
    
    new_cells.append(cell)

if main_loop_found:
    nb["cells"] = new_cells
    with open(notebook_path, "w") as f:
        json.dump(nb, f, indent=1)
    print("✅ Notebook successfully updated.")
else:
    print("❌ Could not find the Main Loop cell to replace.")
