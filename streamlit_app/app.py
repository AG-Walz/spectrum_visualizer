"""MSpecViz Streamlit entry: navigation only — pages live in streamlit_app/views/; 
invoke via streamlit run streamlit_app/app.py from the project root.
"""
import sys
import json
from pathlib import Path
import streamlit as st

# Path setup: put project root on sys.path so import streamlit_app.* and
# import core.* resolve when Streamlit launches this file directly.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# Settings: load app-wide settings.json once per Streamlit session
if "settings" not in st.session_state:
    settings_path = Path(__file__).parent / "config" / "settings.json"
    with open(settings_path, "r") as f:
        st.session_state.settings = json.load(f)

if __name__ == '__main__':
    # Navigation: register the four tabs (Setup, Search & Load, Explore, Export)
    views = Path(__file__).parent / "views"
    pages = {
        str(st.session_state.settings["app_name"]): [
            st.Page(views / "setup_page.py", title="Setup", icon="⚙️"),
            st.Page(views / "search_load" / "search_load_page.py", title="Search & Load", icon="🔍", url_path="search-load"),
            st.Page(views / "explore" / "explore_page.py", title="Explore", icon="🔬", url_path="explore"),
            st.Page(views / "export" / "export_page.py", title="Export", icon="📦", url_path="export"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
