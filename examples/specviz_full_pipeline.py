
"""
SpecViz Full Pipeline - Debugging/Reproduction Script
Converted from examples/specviz_full_pipeline.ipynb
"""

from pathlib import Path
import sys

# from IPython.display import display
def display(obj):
    print(obj)

import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore')

# -------------------------------------------------------------------------
# Setup and Imports
# -------------------------------------------------------------------------

# Add SpecViz paths for imports
cwd = Path.cwd().resolve()
specviz_dir = None

# Try to find specviz directory
if (cwd / "specviz").exists():
    specviz_dir = cwd / "specviz"
elif cwd.name == "specviz":
    specviz_dir = cwd
elif cwd.name == "examples" and cwd.parent.name == "specviz":
    specviz_dir = cwd.parent
else:
    for parent in [cwd] + list(cwd.parents):
        if (parent / "specviz").exists():
            specviz_dir = parent / "specviz"
            break

if specviz_dir is None:
    raise RuntimeError("Could not locate specviz directory. Run from repo root or specviz/examples.")

sys.path.insert(0, str(specviz_dir / "utils"))
sys.path.insert(0, str(specviz_dir / "core"))
sys.path.insert(0, str(specviz_dir / "workspaces"))

print(f"✅ SpecViz directory: {specviz_dir}")

# Import core modules
from shared_imports import pd, oms, Peptidoform

# CRITICAL: Mock Streamlit cache decorators BEFORE importing any modules that use them
from unittest.mock import MagicMock

# Create mock Streamlit module if it doesn't exist or override it
if 'streamlit' not in sys.modules:
    sys.modules['streamlit'] = MagicMock()

# Create no-op cache decorators that just return the function unchanged
def noop_decorator(func=None, **kwargs):
    """No-op decorator for Streamlit cache functions in notebook context."""
    if func is None:
        # Called with arguments: @cache_data(ttl=3600)
        return lambda f: f
    else:
        # Called without arguments: @cache_data
        return func

# Override Streamlit cache decorators with no-ops
sys.modules['streamlit'].cache_data = noop_decorator
sys.modules['streamlit'].cache_resource = noop_decorator

# Now import spectrum_analysis - cache decorators will be no-ops
import spectrum_analysis

# Reload the module to pick up the bug fix for bytes decoding
import importlib
importlib.reload(spectrum_analysis)

# IMPORTANT: Patch the cached loader to avoid pickling pyOpenMS objects in notebooks
def _load_idxml_peptides_no_cache(idxml_path):
    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)
    return protein_ids, peptide_ids

spectrum_analysis._load_idxml_peptides = _load_idxml_peptides_no_cache

# Import the EXACT same pipeline functions used by the app
from spectrum_analysis import (
    filter_idxml_by_selected_peptides,  # PHASE 2: Filter idXMLs by peptide list
    load_only_needed_spectra,           # PHASE 3: Load only needed spectra (memory-efficient)
)

# Import the EXACT same plotting functions used by the app
try:
    from spectrum_engine import (
        create_matplotlib_mirror_plot,      # Create mirror plot figure
        export_matplotlib_mirror_plot,      # Export plot to file
        convert_spectrum_to_df,             # Convert spectrum to DataFrame
        generate_theoretical_spectrum,      # Generate theoretical spectrum
        align_spectra,                      # Align experimental and theoretical spectra
        annotate_peaks,                     # Annotate matched peaks
        extract_annotations,                # Extract annotations for plotting
        EXP_COLORS,                         # Color scheme for experimental spectra
        PRED_COLORS,                        # Color scheme for predicted spectra
    )
except ImportError:
    # Handle if create_matplotlib_mirror_plot is missing by monkey patching or warning
    print("Warning: create_matplotlib_mirror_plot not found in spectrum_engine?")
    from spectrum_engine import *

print("✅ Modules loaded successfully (Streamlit cache decorators disabled for notebook)")


# -------------------------------------------------------------------------
# Helper Functions
# -------------------------------------------------------------------------

def collect_mzml_and_idxml(input_folders):
    mzml_paths = []
    idxml_paths = []
    
    for folder in input_folders:
        folder_path = Path(folder)
        if not folder_path.exists():
            raise FileNotFoundError(f"Folder not found: {folder}")
        
        for path in folder_path.rglob("*"):
            if not path.is_file():
                continue
            
            suffix = path.suffix.lower()
            if suffix == ".mzml":
                mzml_paths.append(str(path))
            elif suffix in {".idxml", ".mzid", ".mzidentml"}:
                idxml_paths.append(str(path))
    
    return sorted(set(mzml_paths)), sorted(set(idxml_paths))

def resolve_reference_idxml(reference_idxml, idxml_paths):
    reference_path = Path(reference_idxml)
    if reference_path.exists():
        resolved = str(reference_path)
        if resolved not in idxml_paths:
            raise ValueError(f"Reference idXML is not in input folders: {reference_idxml}")
        return resolved
    matches = [p for p in idxml_paths if Path(p).name == reference_idxml]
    if not matches:
        # Try fuzzy match
        matches = [p for p in idxml_paths if reference_idxml in Path(p).name]
    if not matches:
        raise ValueError(f"Reference idXML not found by name: {reference_idxml}")
    if len(matches) > 1:
        # Prefer exact name match if multiple fuzzy
        exact = [p for p in matches if Path(p).name == reference_idxml]
        if len(exact) == 1:
            return exact[0]
        raise ValueError(f"Reference idXML name is ambiguous: {reference_idxml}")
    return matches[0]

def build_spectrum_table_from_filtered(idxml_to_mzml_spectra):
    rows = []
    for idxml_name, mzml_map in idxml_to_mzml_spectra.items():
        for mzml_file, spec_id_peptide_pairs in mzml_map.items():
            for item in spec_id_peptide_pairs:
                spec_id = item[0]
                peptidoform = item[1]
                try:
                    unmodified = Peptidoform(peptidoform).sequence
                except Exception:
                    unmodified = ""
                rows.append({
                    "spectrum_id": spec_id,
                    "peptidoform": peptidoform,
                    "peptide_sequence": unmodified,
                    "mzml_file": mzml_file,
                    "idxml_source": idxml_name,
                })
    return pd.DataFrame(rows)

def load_indexed_experiments_for_plotting(mzml_paths, max_workers=4):
    experiments = {}
    print(f"  ⚡ Opening {len(mzml_paths)} files in parallel ({max_workers} workers)...")
    
    def _open_single_exp(path):
        try:
            exp = oms.OnDiscMSExperiment()
            exp.openFile(str(path))
            return Path(path).name, exp, f"  ✅ Opened {Path(path).name} ({exp.getNrSpectra()} spectra)"
        except Exception as e:
            return Path(path).name, None, f"  ❌ Failed to open {Path(path).name}: {e}"

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_open_single_exp, p) for p in mzml_paths]
        for future in as_completed(futures):
            name, exp, msg = future.result()
            print(msg)
            if exp:
                experiments[name] = exp
    return experiments

def build_spectrum_lookups(experiments, max_workers=4):
    lookups = {}
    print(f"  ⚡ Building {len(experiments)} lookups in parallel ({max_workers} workers)...")
    
    def _build_single_lookup(args):
        mzml_name, exp = args
        try:
            first_native_id = exp.getSpectrum(0).getNativeID() if exp.getNrSpectra() > 0 else ""
            if "scan=" in first_native_id:
                pattern = "scan=(?<SCAN>\\d+)"
            elif "index=" in first_native_id:
                pattern = "index=(?<SCAN>\\d+)"
            elif "spectrum=" in first_native_id:
                pattern = "spectrum=(?<SCAN>\\d+)"
            else:
                pattern = "(?<SCAN>\\d+)"
            
            lookup = oms.SpectrumLookup()
            lookup.readSpectra(exp, pattern)
            return mzml_name, lookup, None
        except Exception as exc:
            try:
                native_id_map = {}
                for idx in range(exp.getNrSpectra()):
                    native_id = exp.getSpectrum(idx).getNativeID()
                    if native_id:
                        native_id_map[native_id] = idx
                return mzml_name, native_id_map, None
            except Exception as fallback_exc:
                return mzml_name, None, f"  ⚠️ Failed to build lookup for {mzml_name}: {exc}"

    tasks = list(experiments.items())
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_build_single_lookup, t) for t in tasks]
        for future in as_completed(futures):
            mzml_name, lookup, error_msg = future.result()
            if lookup:
                lookups[Path(mzml_name).stem] = lookup
            if error_msg:
                print(error_msg)
    return lookups

def load_spectrum_from_experiment(experiments, lookups, mzml_file, spectrum_id):
    exp = experiments.get(mzml_file)
    lookup = lookups.get(Path(mzml_file).stem)
    if exp is None or lookup is None:
        return None
    try:
        if isinstance(lookup, dict):
            idx = lookup.get(spectrum_id)
            if idx is None:
                return None
        else:
            idx = lookup.findByNativeID(spectrum_id)
        return exp.getSpectrum(idx)
    except Exception:
        return None

def annotate_spectrum_with_theoretical(exp_spec, peptidoform, charge=None, tol=0.02, included_ions=None):
    if exp_spec is None or exp_spec.size() == 0:
        return pd.DataFrame({"mz": [], "intensity": [], "annotation": []})
    try:
        pf = Peptidoform(str(peptidoform))
    except Exception:
        df = convert_spectrum_to_df(exp_spec)
        if "annotation" not in df.columns:
            df["annotation"] = "unassigned"
        return df[["mz", "intensity", "annotation"]]
    
    try:
        if charge is None or pd.isna(charge):
            max_charge = 2
        else:
            max_charge = int(charge)
    except Exception:
        max_charge = 2
    
    theo_spec = generate_theoretical_spectrum(pf, max_charge=max_charge, included_ions=included_ions)
    theo_df = convert_spectrum_to_df(theo_spec)
    pairs = align_spectra(exp_spec, theo_spec, tol)
    aligned_df = annotate_peaks(exp_spec, theo_df, pairs)
    return extract_annotations(exp_spec, aligned_df)


# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------
if __name__ == "__main__":
    input_folders = [
        "/mnt/volume/nacho/data/HNSCC_lumos_synpeps", 
        "/mnt/volume/nacho/data/HNSCC_timstof_synpeps",
    ]

    peptides = [
        "RTEFNLNQY", "DSELQREGVSHY", "GTDRLLIEY", "DTDTNAQITY", "ASDDFDQASY", 
        "ALLGSAFQL", "ALYDSVILL", "VLQGLEIEL", "GMFGKETYV", "QLLEFNGINL", 
        "ALMDEINFMKM", "KLAEISLGV", "VLVPYEPPQV", "QILDGVHYL", "VYFVAPAKF", 
        "KYIHFTPNF", "NWPSRPYLF", "KYLEKYYNL", "AYNAVSFTW", "KYGVHEAIF", 
        "LTRAGILTF", "SQAFLVSIY", "PVFPTELAF", "ALFAVFQTTF", "TAKSATWTY", 
        "SVYRTVIEY", "VTSSVSSGY", "YVAHERLSY", "LPSAFYIKL", "LPGTAPKLLIY", 
        "LPWAEFESL", "VPSEIDAAF", "SAHDQSLVL", "GEADIMISF", "HEFSSPSHLL", 
        "SETRVLVDI", "KEEPGAFVI", "SEGSPLVVL", "DLRNQLLNYRSAISNHG", 
        "KDTLKKMQKFFGLPQ", "VPRKTWWASKSPDNKPVW", "FPGPNIGGDTHFDSAEPWTVRNEDL", 
        "KDTSIAYFKNKELEQG", "FKRGKLLSYVDAEGNPVG", "GDLELRPSTFLIIDPL", 
        "IQDFIEAEDDLSSFR", "DLRNKILTATVDNANVL", "IRQFTSSSSIKGSSG", 
        "ERPKTSVNLISSLWPTLPSG", "FEQYINNLRRQLDSIVG",
    ]
    
    reference_mzml_files = [
        "211129_SS_HNSCC_synpeps_B35_B40_in_QC_W6-32_20%_DDA#1_400-650mz_msms5.mzML",
        "211206_SS_HNSCC_synpeps_A01_A02_in_QC_W6-32_20%_DDA#1_400-650mz_msms33.mzML",
        "211206_SS_HNSCC_synpeps_A24_B15_in_QC_W6-32_20%_DDA#1_400-650mz_msms35.mzML",
    ]
    
    max_plots = 10
    save_plots = True
    fragment_mass_tolerance = 0.02
    ion_types = ['b', 'y']
    
    plot_output_dir = specviz_dir / "plots" / "full_pipeline"
    analysis_output_dir = specviz_dir / "analysis" / "full_pipeline"
    
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    analysis_output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"   Output directories:\n     - Plots: {plot_output_dir}\n     - Analysis: {analysis_output_dir}")
    
    # -------------------------------------------------------------------------
    # Step 1: Discover Files
    # -------------------------------------------------------------------------
    print("🔍 Discovering files in input folders...")
    mzml_paths, idxml_paths = collect_mzml_and_idxml(input_folders)
    print(f"   Found {len(mzml_paths)} mzML files and {len(idxml_paths)} idXML files")
    
    # Resolve reference mzML files
    reference_mzml_paths = []
    for ref_mzml in reference_mzml_files:
        ref_path = Path(ref_mzml)
        if ref_path.exists():
            reference_mzml_paths.append(str(ref_path))
        else:
            matches = [p for p in mzml_paths if Path(p).name == ref_mzml]
            if matches:
                reference_mzml_paths.append(matches[0])
            else:
                print(f"  ⚠️ Warning: Reference mzML not found: {ref_mzml}")
    
    if not reference_mzml_paths:
        print("  ⚠️ No reference mzML files found. Will attempt to use best-scoring spectra from samples as references.")
    else:
        reference_mzml_names = [Path(p).name for p in reference_mzml_paths]
    
    # -------------------------------------------------------------------------
    # Step 2: Filter idXMLs by Peptide List
    # -------------------------------------------------------------------------
    print("🔬 Filtering idXML files by peptide list...")
    idxml_to_mzml_spectra, filter_messages, filter_empty = filter_idxml_by_selected_peptides(
        idxml_paths,
        mzml_paths,
        peptides,
        progress_callback=None,
    )
    
    if filter_empty:
        print(f"   ⚠️ {len(filter_empty)} files were empty")
    # if filter_messages:
    #     print("   Messages:")
    #     for msg in filter_messages[:5]:
    #         print("     " + msg)
    
    # -------------------------------------------------------------------------
    # Step 3: Peptide Coverage Check
    # -------------------------------------------------------------------------
    print("🔍 Building peptide occurrence table...")
    all_occurrences = []
    for idxml_name, mzml_map in idxml_to_mzml_spectra.items():
        for mzml_file, spec_list in mzml_map.items():
            for item in spec_list:
                spec_id = item[0]
                peptidoform = item[1]
                try:
                    pf = Peptidoform(peptidoform)
                    unmodified = pf.sequence
                    modifications = []
                    for mod in pf.parsed_sequence:
                        if mod[1]:
                            mod_info = str(mod[1])
                            modifications.append(f"{mod[0]}{mod_info}")
                    modifications_str = "; ".join(modifications) if modifications else "None"
                except Exception:
                    unmodified = peptidoform
                    modifications_str = "Modified"
                
                all_occurrences.append({
                    'target_peptide': unmodified,
                    'peptidoform': peptidoform,
                    'modifications': modifications_str,
                    'idxml_file': idxml_name,
                    'mzml_file': mzml_file,
                    'spectrum_id': spec_id,
                })
    
    peptide_occurrence_df = pd.DataFrame(all_occurrences)
    if not peptide_occurrence_df.empty:
        column_order = ['target_peptide', 'peptidoform', 'modifications', 'idxml_file', 'mzml_file', 'spectrum_id']
        peptide_occurrence_df = peptide_occurrence_df[column_order]
        
        output_file = analysis_output_dir / "peptide_occurrences.xlsx"
        peptide_occurrence_df.to_excel(output_file, index=False)
        print(f"✅ Found {len(peptide_occurrence_df)} total spectrum identifications")
    
        # Step 3b: Peptidoform Matrix
        print("🔍 Building peptidoform presence/absence matrix...")
        peptidoform_matrix = peptide_occurrence_df.groupby(['peptidoform', 'idxml_file']).size().unstack(fill_value=0)
        peptidoform_matrix = peptidoform_matrix.applymap(lambda x: '✓' if x > 0 else '✗')
        peptidoform_matrix['Total_Files'] = (peptidoform_matrix == '✓').sum(axis=1)
        
        # Missing peptides
        found_peptides = set(peptide_occurrence_df['target_peptide'].unique())
        missing_peptides = [pep for pep in peptides if pep not in found_peptides]
        
        if missing_peptides:
            idxml_columns = [col for col in peptidoform_matrix.columns if col != 'Total_Files']
            for missing_pep in missing_peptides:
                missing_row = {col: '✗' for col in idxml_columns}
                peptidoform_matrix.loc[missing_pep] = missing_row
                peptidoform_matrix.loc[missing_pep, 'Total_Files'] = 0
    
        peptidoform_matrix = peptidoform_matrix.sort_values(['Total_Files', peptidoform_matrix.index.name], ascending=[False, True])
        
        matrix_output_file = analysis_output_dir / "peptidoform_presence_matrix.xlsx"
        peptidoform_matrix.to_excel(matrix_output_file)
        print(f"   Peptidoforms found: {len(peptidoform_matrix) - len(missing_peptides)}")
    
        spectrum_table = build_spectrum_table_from_filtered(idxml_to_mzml_spectra)
    else:
        print("⚠️ No spectrum identifications found.")
        spectrum_table = pd.DataFrame()
    
    # -------------------------------------------------------------------------
    # Step 4: Load Required Spectra
    # -------------------------------------------------------------------------
    if not spectrum_table.empty:
        print("📂 Loading required spectra from mzML files...")
        needed_mzmls = sorted(set(spectrum_table["mzml_file"].tolist()))
        mzml_paths_filtered = [p for p in mzml_paths if Path(p).name in needed_mzmls]
        
        mapping_df, validation_info, load_messages, load_empty = load_only_needed_spectra(
            idxml_paths,
            mzml_paths_filtered,
            idxml_to_mzml_spectra,
            progress_callback=None,
        )
        
        print(f"✅ Loading complete. Loaded {len(mapping_df)} spectra")
        if not mapping_df.empty:
            mapping_output_file = analysis_output_dir / "spectrum_mapping.xlsx"
            mapping_df.to_excel(mapping_output_file, index=False)
            
        # Summarize errors
        corrupted_files = [m for m in load_messages if "Corrupted" in m or "Failed to open" in m or "Parse Error" in m]
        if corrupted_files:
            print("\n  ❌ CRITICAL ERRORS: The following files could not be loaded:")
            for msg in set(corrupted_files):
                print(f"  {msg.strip()}")
            print("  These files will be excluded from analysis.\n")
    else:
        mapping_df = pd.DataFrame()
    
    # -------------------------------------------------------------------------
    # Step 5: Load mzML Files for Plotting
    # -------------------------------------------------------------------------
    print("🗄️ Opening mzML files with indexed access...")
    experiments = load_indexed_experiments_for_plotting(mzml_paths_filtered)
    lookups = build_spectrum_lookups(experiments)
    print(f"✅ Built {len(lookups)} spectrum lookups")
    
    # -------------------------------------------------------------------------
    # Step 6: Build Reference Lookup
    # -------------------------------------------------------------------------
    print("📚 Building reference lookup...")
    print("📚 Building reference lookup...")
    if not mapping_df.empty:
        reference_rows = mapping_df[mapping_df["mzml_file"].isin(reference_mzml_names)].copy() if reference_mzml_names else pd.DataFrame()
        
        # Auto-Reference Logic
        if reference_rows.empty:
            print("   ⚠️ No external reference spectra found. Selecting best-scoring spectra from samples as internal references...")
            # Heuristic: Sort by score. Assuming higher is better (e.g. Hyperscore, XCorr). 
            # If using e-values/q-values, this might pick the worst. 
            # TODO: robustly determine score direction. For now, assuming standard search engine scores.
            try:
                # Group by peptidoform/sequence and take the one with Max score
                # Note: If score column is string or missing, this might fail slightly or use lexicographical sort
                reference_rows = mapping_df.sort_values("score", ascending=False).groupby("peptide_sequence").head(1)
                print(f"   ℹ️ Selected {len(reference_rows)} internal reference spectra based on max score.")
                
                # Update reference_mzml_names so we know which files to keep loaded/highlight
                reference_mzml_names = reference_rows["mzml_file"].unique().tolist()
            except Exception as e:
                print(f"   ❌ Failed to select auto-references: {e}")
        
    print(f"   Found {len(reference_rows)} reference spectra")
    
    reference_by_unmodified = {}
    for _, row in reference_rows.iterrows():
        try:
            pf = Peptidoform(row["peptide_sequence"])
            unmodified_seq = pf.sequence
        except Exception:
            unmodified_seq = row["peptide_sequence"]
        
        reference_by_unmodified.setdefault(unmodified_seq, []).append(row)
    
    def pick_reference_row(sample_row):
        try:
            pf = Peptidoform(sample_row["peptide_sequence"])
            unmodified_seq = pf.sequence
        except Exception:
            unmodified_seq = sample_row["peptide_sequence"]
        
        candidates = reference_by_unmodified.get(unmodified_seq, [])
        if not candidates:
            return None
        
        sample_charge = sample_row.get("charge")
        if sample_charge is not None and not pd.isna(sample_charge):
            for row in candidates:
                if "charge" in row and row["charge"] == sample_charge:
                    return row
        return candidates[0]
    
    # -------------------------------------------------------------------------
    # Step 7: Generate Mirror Plots
    # -------------------------------------------------------------------------
    print("📊 Generating One-vs-All Mirror Plots...")
    count = 0
    for idx, row in mapping_df.iterrows():
        if max_plots and count >= max_plots:
            break
        
        # Skip reference rows (don't compare reference to itself)
        if row["mzml_file"] in reference_mzml_names:
            continue
            
        ref_row = pick_reference_row(row)
        if ref_row is None:
            continue
        
        # Load spectra
        spec1 = load_spectrum_from_experiment(experiments, lookups, row["mzml_file"], row["spectrum_id"])
        spec2 = load_spectrum_from_experiment(experiments, lookups, ref_row["mzml_file"], ref_row["spectrum_id"])
        
        if spec1 and spec2:
            # Annotate
            df1 = annotate_spectrum_with_theoretical(spec1, row["peptide_sequence"], row["charge"], fragment_mass_tolerance, included_ions=ion_types)
            df2 = annotate_spectrum_with_theoretical(spec2, ref_row["peptide_sequence"], ref_row["charge"], fragment_mass_tolerance, included_ions=ion_types)
            
            # Plot
            title = f"One-vs-All: {row['peptide_sequence']}\nTop: {row['mzml_file']}\nBottom: {ref_row['mzml_file']} (Reference)"
            fig = create_matplotlib_mirror_plot(df1, df2, title=title, r2_metric="both")
            
            if save_plots and fig:
                fname = f"plot_{count}_{row['peptide_sequence']}"
                export_matplotlib_mirror_plot(fig, str(plot_output_dir), fname)
                count += 1
                if count % 5 == 0:
                    print(f"   Generated {count} plots...")
    
        print(f"✅ Generated {count} plots total.")
