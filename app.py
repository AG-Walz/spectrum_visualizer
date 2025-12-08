"""SpecViz App - Main application entry point.

Initializes Streamlit app, preloads packages via shared_imports,
and sets up page navigation. 
"""
import streamlit as st
import sys
from pathlib import Path

# Add utils to path for shared_imports
sys.path.insert(0, str(Path(__file__).parent / "utils"))

# Preload all heavy packages at startup (before user navigates anywhere)
# This ensures all packages (PyOpenMS, psm_utils, etc.) are loaded only once, 
# regardless of workspace selection
with st.spinner("🔄 Loading application modules..."):
    import shared_imports  # noqa: F401

# Now import what we need from shared_imports
from shared_imports import json, os

# Configure Streamlit page
st.set_page_config(
    page_title="SpecViz - MS/MS Spectrum Visualization",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load settings
if "settings" not in st.session_state:
    try:
        settings_path = Path(__file__).parent / "settings.json"
        with open(settings_path, "r") as f:
            st.session_state.settings = json.load(f)
    except FileNotFoundError:
        # Default settings if file doesn't exist
        st.session_state.settings = {
            "app-name": "SpecViz",
            "version": "1.0.0",
            "enable_workspaces": True
        }

# Configure environment
os.environ.setdefault("BOKEH_ALLOW_WS_ORIGIN", "*")

if __name__ == '__main__':
    # Define navigation structure with proper spacing
    app_name = st.session_state.settings.get("app-name", "SpecViz")
    
    pages = {
        f"🏠 {app_name}": [
            st.Page(Path("pages", "quickstart.py"), title="Quickstart", icon="🏁"),
            st.Page(Path("pages", "documentation.py"), title="Documentation", icon="📖"),
        ],
        "🔬 Workspaces": [
            st.Page(Path("workspaces", "spectrum_analysis.py"), title="Spectrum Analysis", icon="📊"),
            st.Page(Path("workspaces", "dataset_filter.py"), title="Filter Large Datasets", icon="🗂️"),
        ],
        "⚙️ Settings": [
            st.Page(Path("workspaces", "settings.py"), title="Configuration", icon="🛠️"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()