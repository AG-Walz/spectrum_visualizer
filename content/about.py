"""
About MSpecViz page.
"""
import streamlit as st

from src.common.common import page_setup

params = page_setup()

st.title("About MSpecViz")

st.markdown("""
**MSpecViz** is a mass spectrometry spectrum visualization tool for proteomics research.

### Features
- Search peptides across multiple idXML files
- Memory-efficient mzML spectrum loading
- Mirror plots comparing query vs reference spectra
- Peak annotation with theoretical fragment ions
- Export to PNG/SVG with organized folder structure

### Usage
1. Select idXML identification files
2. Select corresponding mzML spectrum files
3. Enter target peptides (one per line)
4. Mark reference files using substring filter
5. Configure ion types and m/z tolerance
6. Generate and export mirror plots
""")

st.markdown("---")

# Configuration display
with st.expander("Current Configuration"):
    st.json(st.session_state.get("settings", {}))
