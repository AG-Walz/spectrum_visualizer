import streamlit as st
import json
from pathlib import Path

if "settings" not in st.session_state:
    settings_path = Path(__file__).parent / "settings.json"
    with open(settings_path, "r") as f:
        st.session_state.settings = json.load(f)

if __name__ == '__main__':
    pages = {
        str(st.session_state.settings["app-name"]): [
            st.Page(Path("content", "setup.py"), title="Setup", icon="⚙️"),
            st.Page(Path("content", "search_load.py"), title="Search & Load", icon="🔍"),
            st.Page(Path("content", "explore.py"), title="Explore", icon="🔬"),
            st.Page(Path("content", "export.py"), title="Export", icon="📦"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
