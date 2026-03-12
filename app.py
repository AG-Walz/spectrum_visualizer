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
            st.Page(Path("content", "peptide_mirror.py"), title="Peptide Mirror Analysis", icon="🔬"),
            st.Page(Path("content", "about.py"), title="About", icon="📖"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
