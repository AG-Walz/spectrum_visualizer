import streamlit as st

# Simple setup - no need for complex common.py
st.set_page_config(
    page_title="immunuXplore - Quickstart",
    initial_sidebar_state="expanded"
)

# Load and display README or custom quickstart content
try:
    with open("README.md", "r") as f:
        readme = f.read()
    st.markdown(readme, unsafe_allow_html=True)
except FileNotFoundError:
    st.title("🔬 SpecViz - Quickstart")
    st.markdown("Welcome to **SpecViz**, an MS/MS spectrum visualization tool")

    st.markdown("---")  # Separator
    
    st.markdown("""
    ## 🏁 Quick Start Guide
    
    #### 1. Navigate to Workspaces
    Once you have read the documentation, click on **"Workspaces"** in the sidebar to access the main spectrum analysis interface. If you exit the workspace, your progress will be resetted, so make sure to read the documentation properly or write the steps down before proceeding.

    #### 2. Upload Your Data
    - **Select mzML files**: Upload all MS/MS spectrum files (.mzML) corresponding to the idXML file using the browse files button. Select multiple at once by pressing shift while clicking on the selection. Then click *Open*.
    - **Select idXML file**: Upload peptide identifications file (.idXML) using the browse files button. Only one file is allowed. Once file has been browsed, click then click *Select*.
    - **Process files**: Once two green checkmarks ✅ appear, click the *Process Files* button to check for any missing or unmatching files.

    #### 3. Select Analysis Parameters
    - **Select a Peptide:** Select the peptide you want to analyze from the dropdown menu.
    - **Choose an AI model:**  For spectrum prediction choose from the dropdown menu between MS2PiP, Prosit, etc. (recommended MS2PIP TimsTOF 2024)
    - **Set Parameters:** Type in the Maximum Fragment Charge (default: 2) and the Mass to Charge Tolerance (default: 0.5)
                
    #### 4. Start analysis
    - Click on **"Start Analysis"** to begin the processing of your data.
    - If a peptide with only one corresponding spectrum was selected, then a mirror plot with the corresponding predicted spectrum and the peaks coloured by OpenMS's TheoreticalSpectrumGenerator will be displayed. If multiple spectra are available, an additional comparison view will be shown.

    """)