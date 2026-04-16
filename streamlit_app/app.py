import sys
import json
from pathlib import Path

# Ensure the project root is on sys.path so that `import streamlit_app.*`
# (and `import mspecviz.*`) resolves when Streamlit runs this file directly.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import streamlit as st

if "settings" not in st.session_state:
    settings_path = Path(__file__).parent / "config" / "settings.json"
    with open(settings_path, "r") as f:
        st.session_state.settings = json.load(f)

if __name__ == '__main__':
    views = Path(__file__).parent / "views"
    pages = {
        str(st.session_state.settings["app_name"]): [
            st.Page(views / "setup.py", title="Setup", icon="⚙️"),
            st.Page(views / "search_load" / "page.py", title="Search & Load", icon="🔍"),
            st.Page(views / "explore" / "page.py", title="Explore", icon="🔬"),
            st.Page(views / "export" / "page.py", title="Export", icon="📦"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
