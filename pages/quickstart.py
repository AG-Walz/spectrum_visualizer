"""Quickstart guide - Getting started for new users.

Step-by-step workflow overview with instructions. Static content page.
Accessed via: App sidebar → Home → Quickstart
"""
import streamlit as st
from pathlib import Path

# Simple setup
st.set_page_config(
    page_title="SpecViz - Quickstart",
    initial_sidebar_state="expanded"
)

# Load markdown documentation (single source of truth)
docs_path = Path(__file__).parent.parent / "docs" / "QUICKSTART.md"

try:
    with open(docs_path, "r", encoding="utf-8") as f:
        quickstart_content = f.read()
    st.markdown(quickstart_content, unsafe_allow_html=True)
except FileNotFoundError:
    st.error(f"❌ Documentation file not found: {docs_path}")
    st.info("Please ensure docs/QUICKSTART.md exists in the project root.")
