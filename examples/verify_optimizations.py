import sys
import time
from pathlib import Path
import os

# Setup paths (assuming running from examples dir or with specviz in path)
cwd = Path.cwd()
specviz_dir = None
if (cwd / "specviz").exists():
    specviz_dir = cwd / "specviz"
elif cwd.name == "examples" and cwd.parent.name == "specviz":
    specviz_dir = cwd.parent
else:
    # Try absolute path
    specviz_dir = Path("/mnt/volume/nacho/specviz")

if not specviz_dir.exists():
    print("SpecViz dir not found")
    sys.exit(1)

sys.path.insert(0, str(specviz_dir / "utils"))
sys.path.insert(0, str(specviz_dir / "core"))
sys.path.insert(0, str(specviz_dir / "workspaces"))
sys.path.insert(0, str(specviz_dir / "examples")) # To import specviz_full_pipeline

# Mock streamlit before imports
from unittest.mock import MagicMock
if 'streamlit' not in sys.modules:
    sys.modules['streamlit'] = MagicMock()
def noop_decorator(func=None, **kwargs):
    if func is None: return lambda f: f
    return func
sys.modules['streamlit'].cache_data = noop_decorator
sys.modules['streamlit'].cache_resource = noop_decorator

from shared_imports import oms # Check import works

# Import optimized functions
import spectrum_analysis
# Patch loader as per pipeline
def _load_idxml_peptides_no_cache(idxml_path):
    protein_ids = []
    peptide_ids = []
    oms.IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)
    return protein_ids, peptide_ids
spectrum_analysis._load_idxml_peptides = _load_idxml_peptides_no_cache

from spectrum_analysis import filter_idxml_by_selected_peptides, load_only_needed_spectra

# Import pipeline helpers
# Note: since we run as script, we can import from the file if we are in same dir, or exec it.
# Easier to just verify spectrum_analysis first, then copy the logic for testing helpers if needed.
# But I want to verify specviz_full_pipeline.py explicitly.
try:
    from specviz_full_pipeline import load_indexed_experiments_for_plotting, build_spectrum_lookups
except ImportError:
    # If define in same dir it works
    sys.path.append(str(specviz_dir / "examples"))
    from specviz_full_pipeline import load_indexed_experiments_for_plotting, build_spectrum_lookups

# Test Data
data_dir = Path("/mnt/volume/nacho/data/HNSCC_lumos_synpeps")
idxml_path = list(data_dir.glob("*_filtered.idXML"))[0]
mzml_paths = [str(p) for p in data_dir.glob("*.mzML")]

# Extract some peptides from the idXML to ensure we get matches
print("\n--- Extracting valid peptides for testing ---")
prot_ids, pep_ids = _load_idxml_peptides_no_cache(idxml_path)
valid_peptides = set()
for p in pep_ids[:100]: # Check first 100
    hits = p.getHits()
    if hits:
        valid_peptides.add(hits[0].getSequence().toString())
        valid_peptides.add(hits[0].getSequence().toUnmodifiedString())

peptides = list(valid_peptides)[:20]
print(f"Selected {len(peptides)} valid peptides for testing")

print(f"Running verification with:\n ID: {idxml_path}\n mzMLs: {len(mzml_paths)}")

# 1. Filter Phase
print("\n--- Testing Filter Phase ---")
start = time.time()
idxml_to_mzml_spectra, msgs, empty = filter_idxml_by_selected_peptides([str(idxml_path)], mzml_paths, peptides)
end = time.time()
print(f"Filter phase took {end-start:.4f}s")
count = sum(len(v) for m in idxml_to_mzml_spectra.values() for v in m.values())
print(f"Found {count} spectra matches")
if count > 0:
    first_key = list(idxml_to_mzml_spectra.keys())[0]
    first_mzml = list(idxml_to_mzml_spectra[first_key].keys())[0]
    print(f"Sample item: {idxml_to_mzml_spectra[first_key][first_mzml][0]}")

# 2. Load Phase
print("\n--- Testing Load Phase ---")
start = time.time()
mapping_df, val_info, msgs, empty = load_only_needed_spectra([str(idxml_path)], mzml_paths, idxml_to_mzml_spectra)
end = time.time()
print(f"Load phase took {end-start:.4f}s")
print(f"Loaded {len(mapping_df)} spectra")

# 3. Pipeline Helpers
if not mapping_df.empty:
    print("\n--- Testing Pipeline Helpers (Experiments & Lookups) ---")
    needed_mzmls = mapping_df["mzml_file"].unique()
    load_paths = [p for p in mzml_paths if Path(p).name in needed_mzmls]
    
    start = time.time()
    experiments = load_indexed_experiments_for_plotting(load_paths)
    print(f"Opened {len(experiments)} experiments in {time.time()-start:.4f}s")
    
    start = time.time()
    lookups = build_spectrum_lookups(experiments)
    print(f"Built {len(lookups)} lookups in {time.time()-start:.4f}s")
