"""Comprehensive documentation - User guide and API reference.

Covers file formats, AI models, workflows, troubleshooting, and technical details.
Accessed via: App sidebar → Home → Documentation

All documentation content is loaded from docs/DOCUMENTATION.md
"""
import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="SpecViz - Documentation",
    initial_sidebar_state="expanded"
)

# Load markdown documentation (single source of truth)
docs_path = Path(__file__).parent.parent / "docs" / "DOCUMENTATION.md"

try:
    with open(docs_path, "r", encoding="utf-8") as f:
        documentation_content = f.read()
    st.markdown(documentation_content, unsafe_allow_html=True)
except FileNotFoundError:
    st.error(f"❌ Documentation file not found: {docs_path}")
