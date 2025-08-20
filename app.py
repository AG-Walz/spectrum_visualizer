import streamlit as st
from pathlib import Path
import json
import os

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
        with open("settings.json", "r") as f:
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
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="🏁"),
            st.Page(Path("content", "documentation.py"), title="Documentation", icon="📖"),
        ],
        "🔬 Workspaces": [
            st.Page(Path("src", "spectrum_workspace.py"), title="Spectrum Analysis", icon="📊"),
        ],
        "⚙️ Settings": [
            st.Page(Path("src", "settings.py"), title="Configuration", icon="🛠️"),
        ],
        "ℹ️ About": [
            st.Page(Path("content", "version.py"), title="Version & Info", icon="📜"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()