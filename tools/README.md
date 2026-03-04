# MSpecViz Tools

## convert_mzml_to_parquet.py

Standalone conversion script that extracts matched spectra from mzML files (based on idXML search results) and writes them to a compact Parquet file. This dramatically speeds up repeated analysis — loading from Parquet is ~100x faster than re-parsing mzML.

### Requirements

- Python 3.11+
- pyopenms
- pyarrow
- numpy

These are all included in the MSpecViz conda environment:
```bash
conda activate mspecviz
```

The script also imports from the `src/` modules in the parent directory (peptide search and spectrum loading). When running from the project root, this is handled automatically.

### Quick Start

```bash
# From the MSpecViz project root:
conda run -n mspecviz python tools/convert_mzml_to_parquet.py \
  --idxml /path/to/file1.idXML /path/to/file2.idXML \
  --mzml-dir /path/to/mzml_files/ \
  --peptides "PEPTIDEK,ANOTHERSEQ" \
  --output /path/to/output.parquet
```

### CLI Usage

```
usage: convert_mzml_to_parquet.py [-h] --idxml IDXML [IDXML ...]
                                   --mzml-dir MZML_DIR
                                   --peptides PEPTIDES
                                   --output OUTPUT
                                   [--reference-substring REFERENCE_SUBSTRING]
                                   [--verbose]

Arguments:
  --idxml           Path(s) to idXML identification files
  --mzml-dir        Directory containing mzML spectrum files
  --peptides        Comma-separated target peptide sequences (base sequences)
  --output          Output Parquet file path
  --reference-substring  Substring to identify reference mzML files (e.g., "Synpep")
  -v, --verbose     Enable debug logging
```

### Python API

The script can also be imported and used programmatically:

```python
from tools.convert_mzml_to_parquet import convert_mzml_to_parquet

output_path = convert_mzml_to_parquet(
    idxml_paths=["data/search1.idXML", "data/search2.idXML"],
    mzml_dir="data/mzml/",
    target_peptides=["PEPTIDEK", "ANOTHERSEQ"],
    output_path="output/spectra.parquet",
    reference_substring="Synpep",          # Optional
    progress_callback=print,               # Optional
    error_callback=lambda e: print(f"ERR: {e}"),  # Optional
)
```

### What It Does

1. **Searches idXML files** for target peptides using PyOpenMS
2. **Resolves mzML filenames** via `id_merge_index` metadata (correct provenance mapping)
3. **Opens each mzML** using `OnDiscMSExperiment` (memory-efficient, no full load)
4. **Extracts only matching spectra** using optimized targeted index lookup
5. **Validates** each spectrum (MS level check, provenance verification)
6. **Writes Parquet** with all data needed for downstream analysis

### Output Parquet Schema

Each row represents one spectrum:

| Column | Type | Description |
|--------|------|-------------|
| source_mzml | string | Original mzML filename (provenance) |
| spectrum_id | string | Native ID (e.g., "index=123") |
| peptidoform | string | Full sequence with modifications |
| base_sequence | string | Unmodified peptide sequence |
| variant_key | string | Grouping key (sequence\|modification) |
| charge | int32 | Precursor charge state |
| precursor_mz | float64 | Precursor m/z |
| rt | float64 | Retention time (seconds) |
| is_reference | bool | True if from a reference file |
| mz | list\<float64\> | m/z array (variable length) |
| intensity | list\<float64\> | Intensity array (variable length) |
| score | float64 | PSM score |
| idxml_file | string | Source idXML filename |

The file also contains Parquet metadata with provenance info (conversion timestamp, source files, warning counts).

### Using the Output

**In MSpecViz** — Select the .parquet file in the Files tab and go straight to Process. Skips search and mzML parsing entirely.

**In Python:**
```python
import pandas as pd
df = pd.read_parquet("output.parquet")
print(df[["base_sequence", "source_mzml", "is_reference"]].head())

# Access spectrum data
for _, row in df.iterrows():
    mz = row["mz"]        # list of floats
    intensity = row["intensity"]  # list of floats
```

**In R:**
```r
library(arrow)
df <- read_parquet("output.parquet")
```

### Example with Data

```bash
# Convert a dataset to Parquet
conda run -n mspecviz python tools/convert_mzml_to_parquet.py \
  --idxml /data/my_project/identifications/*.idXML \
  --mzml-dir /data/my_project/spectra/ \
  --peptides "LGGNEQVTR,GAGSSEPVTGLDAK" \
  --output exports/my_spectra.parquet \
  --reference-substring "Reference" \
  -v

# Verify output
conda run -n mspecviz python -c "
import pyarrow.parquet as pq
t = pq.read_table('exports/my_spectra.parquet')
print(t.schema)
print(f'Rows: {len(t)}')
print(t.to_pandas()[['source_mzml','spectrum_id','base_sequence','is_reference']].head(20))
"
```
