
import sys
import os
from pyopenms import IdXMLFile, ProteinIdentification, PeptideIdentification

# Path to an idXML file (adjust if necessary)
# Based on notebook output: /mnt/volume/nacho/data/ProCa_Anna/idxml_classI/ProC_synpeps_a_pout_filtered.idXML
idxml_path = "/mnt/volume/nacho/data/ProCa_Anna/idxml_classI/ProC_synpeps_a_pout_filtered.idXML"

if not os.path.exists(idxml_path):
    print(f"File not found: {idxml_path}")
    # Try searching for one
    import glob
    files = glob.glob("/mnt/volume/nacho/data/ProCa_Anna/idxml_classI/*.idXML")
    if files:
        idxml_path = files[0]
        print(f"Using found file: {idxml_path}")
    else:
        sys.exit(1)

print(f"Inspecting: {idxml_path}")

protein_ids = []
peptide_ids = []
IdXMLFile().load(idxml_path, protein_ids, peptide_ids)

print(f"ProteinIdentification objects: {len(protein_ids)}")
for i, pid in enumerate(protein_ids):
    print(f"  [{i}] Identifier: {pid.getIdentifier()}")
    # Check for metadata
    keys = ["CometAdapter:1:in", "spectra_data", "spectrum_reference", "mzML_file"]
    for key in keys:
        if pid.metaValueExists(key):
             print(f"      {key}: {pid.getMetaValue(key)}")

print(f"\nPeptideIdentification objects: {len(peptide_ids)}")
print("First 5 peptides:")
for i, pep in enumerate(peptide_ids[:5]):
    print(f"  [{i}] Identifier: {pep.getIdentifier()}")
    print(f"      Spectrum Ref: {pep.getMetaValue('spectrum_reference') if pep.metaValueExists('spectrum_reference') else 'N/A'}")
    print(f"      File Origin: {pep.getMetaValue('file_origin') if pep.metaValueExists('file_origin') else 'N/A'}")
    print(f"      mzML_file: {pep.getMetaValue('mzML_file') if pep.metaValueExists('mzML_file') else 'N/A'}")
    print(f"      spectra_data: {pep.getMetaValue('spectra_data') if pep.metaValueExists('spectra_data') else 'N/A'}")

# Check if identifiers match
pid_identifiers = set(p.getIdentifier() for p in protein_ids)
pep_identifiers = set(p.getIdentifier() for p in peptide_ids)

print(f"\nUnique ProteinID Identifiers: {pid_identifiers}")
print(f"Unique PeptideID Identifiers: {pep_identifiers}")
