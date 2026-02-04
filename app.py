"""
MSpecViz - MS/MS Spectrum Visualization

Main Streamlit entry point.
"""
import streamlit as st
import json
from pathlib import Path


# Page configuration - must be first Streamlit call
st.set_page_config(
    page_title="MSpecViz",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_settings():
    """Load application settings from settings.json."""
    if 'settings' not in st.session_state:
        settings_path = Path(__file__).parent / "settings.json"
        if settings_path.exists():
            with open(settings_path) as f:
                st.session_state.settings = json.load(f)
        else:
            st.session_state.settings = {}
    return st.session_state.settings


def main():
    """Main application entry point."""
    settings = load_settings()

    # Sidebar navigation
    st.sidebar.title(f"🔬 {settings.get('app_name', 'MSpecViz')}")
    st.sidebar.caption(f"v{settings.get('version', '0.1.0')}")

    # Workspace selection
    workspace = st.sidebar.radio(
        "Workspace",
        options=["Peptide Mirror Analysis", "About"],
        label_visibility="collapsed"
    )

    st.sidebar.markdown("---")

    # Main content area
    if workspace == "Peptide Mirror Analysis":
        from workspaces import peptide_mirror
        peptide_mirror.render()
    elif workspace == "About":
        render_about(settings)


def render_about(settings):
    """Render the About page."""
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
        st.json(settings)


if __name__ == "__main__":
    main()
