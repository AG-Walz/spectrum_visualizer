import streamlit as st
import json
from pathlib import Path
import datetime

def page_setup():
    """Common page setup function."""
    if "settings" not in st.session_state:
        try:
            with open("settings.json", "r") as f:
                st.session_state.settings = json.load(f)
        except FileNotFoundError:
            st.session_state.settings = {
                "app-name": "SpecViz",
                "version": "0.0.1"
            }

page_setup()

# App information
app_name = st.session_state.settings.get("app-name", "SpecViz")
version = st.session_state.settings.get("version", "0.0.1")

st.title(f"ℹ️ About {app_name}")

st.markdown("---")  # Separator

# Main info section
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown(f"""
    ## 🔬 {app_name} v{version}
    
    **MS/MS Spectrum Visualization Platform**
    
    {app_name} is a streamlit application for visualizing mass spectrometry data with 
    integrated koina predictions. 
    
    ### ✨ Key Features
    
    - **Interactive Visualization**: Mirror plots comparing experimental vs predicted spectra
    - **File Format Support**: mzML, idXML, mzIdentML compatibility
    - **Peptide Modifications**: Full support with proper mass shift formatting
    - **Modern Interface**: Clean, responsive design with comprehensive help
    - **Robust Processing**: Error handling, validation, and multiple fallback methods
    
    ### 🎯 Use Cases
    
    - **Proteomics Research**: Validate peptide identifications
    - **Method Development**: Compare Koina model performance
    - **Quality Control**: Assess spectrum quality and annotations
    - **Publication**: Generate spectrum figures & stuff
    """)

with col2:
    # Logo/icon area
    st.markdown("""
    <div style="text-align: center; padding: 20px;">
        <div style="font-size: 80px;">🔬</div>
        <h3>SpecViz</h3>
        <p><em>Spectrum Visualization</em></p>
    </div>
    """, unsafe_allow_html=True)
    
    # Quick stats
    st.markdown("### 📊 Quick Stats")
    st.metric("Version", version)
    st.metric("Release", "July 2025")
    st.metric("Platform", "Streamlit")

# Technical Information
st.markdown("---")  # Separator
st.header("🔧 Technical Information")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Architecture")
    st.markdown("""
    - **Frontend**: Streamlit web framework
    - **Backend**: OpenMS for MS data processing  
    - **AI Integration**: Koina API for model access
    - **Visualization**: Bokeh for interactive plots
    - **Data Processing**: Pandas and NumPy
    """)
    
    st.subheader("Supported Formats")
    st.markdown("""
    **Input Files:**
    - mzML (mass spectrometry data)
    - idXML (OpenMS identifications)
    - mzIdentML (PSI standard)
    
    **AI Models:**
    - MS2PiP (timsTOF2024, HCD2021)
    - Prosit (2019, 2020 intensity)
    - AlphaPept (generic MS2)
    """)

with col2:
    st.subheader("Dependencies")
    
    # Try to get package versions
    try:
        import streamlit as st_version
        import pandas as pd_version
        import numpy as np_version
        
        dependencies = {
            "Streamlit": st_version.__version__,
            "Pandas": pd_version.__version__,
            "NumPy": np_version.__version__,
        }
        
        try:
            import bokeh
            dependencies["Bokeh"] = bokeh.__version__
        except ImportError:
            dependencies["Bokeh"] = "Not installed"
        
        try:
            import requests
            dependencies["Requests"] = requests.__version__
        except ImportError:
            dependencies["Requests"] = "Not installed"
            
        for name, ver in dependencies.items():
            st.code(f"{name}: {ver}")
            
    except Exception:
        st.info("Version information not available")
    
    st.subheader("System Requirements")
    st.markdown("""
    - Python 3.8+
    - 4GB+ RAM recommended
    - Modern web browser
    - Internet connection (for AI models)
    """)
st.markdown("---")  # Separator
# Development Information
st.header("👥 Development & Support")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Development Team")
    st.markdown("""
    **Spectrum Visualization Team**
    
    - Modern UI/UX design
    - Robust data processing
    - AI model integration
    - Comprehensive testing
    
    *Built with ❤️ for the proteomics community*
    """)
    
    st.subheader("Acknowledgments")
    st.markdown("""
    - **OpenMS**: Mass spectrometry data processing
    - **Koina**: AI model serving platform
    - **Streamlit**: Web application framework
    - **Proteomics Community**: Feedback and testing
    """)

with col2:
    st.subheader("Support & Resources")
    
    # GitHub link (if available)
    github_user = st.session_state.settings.get("github-user", "spectrum-visualization-team")
    repo_name = st.session_state.settings.get("repository-name", "specviz")
    
    st.markdown(f"""
    - 📚 [Documentation](#{app_name.lower()}-documentation)
    - 🐛 [Report Issues](https://github.com/{github_user}/{repo_name}/issues)
    - 💡 [Feature Requests](https://github.com/{github_user}/{repo_name}/discussions)
    - 📧 [Contact Support](mailto:support@{app_name.lower()}.org)
    """)
    
    st.subheader("Contributing")
    st.markdown("""
    We welcome contributions!
    
    - Submit bug reports
    - Suggest new features  
    - Improve documentation
    - Add new AI models
    """)

st.markdown("---")  # Separator
# License and Legal
st.header("📜 License & Legal")

col1, col2 = st.columns(2)

with col1:
    st.subheader("License")
    st.markdown("""
    **MIT License**
    
    Copyright (c) 2025 Spectrum Visualization Team
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software.
    """)

with col2:
    st.subheader("Privacy & Data")
    st.markdown("""
    **Data Processing:**
    - Files processed locally/temporarily
    - No data stored permanently
    - AI predictions via Koina API
    
    **Privacy:**
    - No personal data collection
    - Optional anonymous analytics
    - Configurable in Settings
    """)

# Changelog
st.header("📋 Changelog")

with st.expander("Version History"):
    st.markdown(f"""
    ### v{version} (Current)
    - ✨ Complete workspace structure implementation
    - 🔧 Settings and configuration management
    - 📖 Comprehensive documentation
    - 🎨 Modern UI with immunoexplore-style navigation
    - 🤖 Multiple AI model support
    - 📊 Resource monitoring and performance metrics
    - 🐛 Robust error handling and debugging
    
    ### v0.9.0 (Previous)
    - 🔬 Core spectrum visualization functionality
    - 📁 File upload and validation
    - 🎯 Basic AI model integration
    - 📈 Bokeh plotting with fallbacks
    
    ### v0.1.0 (Initial)
    - 🏗️ Basic Streamlit application structure
    - 📊 Simple spectrum display
    - 🔍 File format support
    """)

# Footer
st.markdown("---")
st.markdown(f"""
<div style="text-align: center; color: #666; margin-top: 40px;">
    <p>{app_name} v{version} • Built with Streamlit • {datetime.datetime.now().year}</p>
    <p><em>Advancing proteomics through AI-powered visualization</em></p>
</div>
""", unsafe_allow_html=True)