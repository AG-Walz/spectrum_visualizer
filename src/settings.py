import streamlit as st
import json
from pathlib import Path

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

st.title("⚙️ Settings & Configuration")
st.markdown("Configure SpecViz application settings and preferences.")

st.markdown("---")  # Separator

# Load current settings
settings = st.session_state.settings

# Application Settings
st.header("🔧 Application Settings")

col1, col2 = st.columns(2)

with col1:
    st.subheader("General")
    
    app_name = st.text_input(
        "Application Name",
        value=settings.get("app-name", "SpecViz"),
        help="Display name for the application"
    )
    
    version = st.text_input(
        "Version",
        value=settings.get("version", "0.0.1"),
        help="Application version number"
    )
    
    enable_workspaces = st.checkbox(
        "Enable Workspaces",
        value=settings.get("enable_workspaces", True),
        help="Enable workspace functionality"
    )

with col2:
    st.subheader("User Interface")
    
    theme = st.selectbox(
        "Theme",
        options=["light", "dark", "auto"],
        index=["light", "dark", "auto"].index(settings.get("ui", {}).get("theme", "light")),
        help="Application color theme"
    )
    
    sidebar_expanded = st.checkbox(
        "Sidebar Expanded by Default",
        value=settings.get("ui", {}).get("sidebar_expanded", True),
        help="Default sidebar state"
    )
    
    debug_mode_default = st.checkbox(
        "Debug Mode Default",
        value=settings.get("ui", {}).get("debug_mode", False),
        help="Enable debug mode by default"
    )

# Model Settings
st.markdown("---")  # Separator
st.header("🤖 AI Model Configuration")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Default Model")
    
    available_models = [
        "ms2pip_timsTOF2024",
        "ms2pip_HCD2021", 
        "Prosit_2019_intensity",
        "Prosit_2020_intensity",
        "AlphaPept_ms2_generic"
    ]
    
    current_default = settings.get("models", {}).get("default", "ms2pip_timsTOF2024")
    default_model = st.selectbox(
        "Default AI Model",
        options=available_models,
        index=available_models.index(current_default) if current_default in available_models else 0,
        help="Default model for spectrum prediction"
    )

with col2:
    st.subheader("Available Models")
    
    # Multi-select for available models
    current_available = settings.get("models", {}).get("available", available_models)
    selected_models = st.multiselect(
        "Enabled Models",
        options=available_models,
        default=current_available,
        help="Models available for selection"
    )

# Performance Settings
st.markdown("---")  # Separator
st.header("⚡ Performance Settings")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Resource Monitoring")
    
    enable_monitoring = st.checkbox(
        "Enable Resource Monitoring",
        value=True,
        help="Monitor CPU and memory usage"
    )
    
    cache_size = st.slider(
        "Cache Size (MB)",
        min_value=100,
        max_value=2000,
        value=settings.get("performance", {}).get("cache_size_mb", 500),
        step=100,
        help="Maximum cache size for file operations"
    )

with col2:
    st.subheader("File Processing")
    
    max_file_size = st.slider(
        "Max File Size (GB)",
        min_value=1,
        max_value=10,
        value=settings.get("performance", {}).get("max_file_size_gb", 2),
        help="Maximum allowed file upload size (requires app restart to apply)"
    )
    
    chunk_size = st.slider(
        "Processing Chunk Size",
        min_value=1000,
        max_value=10000,
        value=settings.get("performance", {}).get("chunk_size", 5000),
        step=1000,
        help="Number of spectra processed per chunk"
    )

# File upload information
st.info(
    "💡 **File Upload Limits**: "
    f"Current maximum upload size is set to {max_file_size} GB. "
    "This setting is also configured in `.streamlit/config.toml`. "
    "Changes require an application restart to take effect."
)

# Analytics Settings
st.markdown("---")  # Separator
st.header("📊 Analytics & Monitoring")

enable_analytics = st.checkbox(
    "Enable Usage Analytics",
    value=settings.get("analytics", {}).get("enabled", False),
    help="Collect anonymous usage statistics"
)

if enable_analytics:
    col1, col2 = st.columns(2)
    
    with col1:
        ga_enabled = st.checkbox(
            "Google Analytics",
            value=settings.get("analytics", {}).get("google-analytics", {}).get("enabled", False)
        )
        
        if ga_enabled:
            ga_tag = st.text_input(
                "GA Tag",
                value=settings.get("analytics", {}).get("google-analytics", {}).get("tag", ""),
                help="Google Analytics tracking ID"
            )
    
    with col2:
        piwik_enabled = st.checkbox(
            "Piwik Pro",
            value=settings.get("analytics", {}).get("piwik-pro", {}).get("enabled", False)
        )
        
        if piwik_enabled:
            piwik_tag = st.text_input(
                "Piwik Tag",
                value=settings.get("analytics", {}).get("piwik-pro", {}).get("tag", ""),
                help="Piwik Pro tracking ID"
            )

# Advanced Settings
st.markdown("---")  # Separator
st.header("🔬 Advanced Configuration")

with st.expander("Developer Settings"):
    test_mode = st.checkbox(
        "Test Mode",
        value=settings.get("test", False),
        help="Enable test mode for development"
    )
    
    workspaces_dir = st.text_input(
        "Workspaces Directory",
        value=settings.get("workspaces_dir", ".."),
        help="Directory for workspace data storage"
    )
    
    online_deployment = st.checkbox(
        "Online Deployment",
        value=settings.get("online_deployment", False),
        help="Enable features for online deployment"
    )

# Save Settings
st.markdown("---")  # Separator
st.header("💾 Save Configuration")

col1, col2, col3 = st.columns([1, 1, 2])

with col1:
    if st.button("💾 Save Settings", type="primary"):
        # Update settings
        new_settings = {
            "app-name": app_name,
            "version": version,
            "enable_workspaces": enable_workspaces,
            "ui": {
                "theme": theme,
                "sidebar_expanded": sidebar_expanded,
                "debug_mode": debug_mode_default
            },
            "models": {
                "default": default_model,
                "available": selected_models
            },
            "analytics": {
                "enabled": enable_analytics,
                "google-analytics": {
                    "enabled": ga_enabled if enable_analytics else False,
                    "tag": ga_tag if enable_analytics and ga_enabled else ""
                },
                "piwik-pro": {
                    "enabled": piwik_enabled if enable_analytics else False,
                    "tag": piwik_tag if enable_analytics and piwik_enabled else ""
                }
            },
            "performance": {
                "cache_size_mb": cache_size,
                "max_file_size_gb": max_file_size,
                "chunk_size": chunk_size,
                "enable_monitoring": enable_monitoring
            },
            "advanced": {
                "test": test_mode,
                "workspaces_dir": workspaces_dir,
                "online_deployment": online_deployment
            }
        }
        
        # Save to file
        try:
            with open("settings.json", "w") as f:
                json.dump(new_settings, f, indent=4)
            
            # Update session state
            st.session_state.settings = new_settings
            
            st.success("✅ Settings saved successfully!")
            st.info("🔄 Restart the application to apply all changes.")
            
        except Exception as e:
            st.error(f"❌ Failed to save settings: {e}")

with col2:
    if st.button("🔄 Reset to Defaults"):
        default_settings = {
            "app-name": "SpecViz",
            "version": "0.0.1",
            "enable_workspaces": True,
            "ui": {
                "theme": "dark",
                "sidebar_expanded": True,
                "debug_mode": False
            },
            "models": {
                "default": "ms2pip_timsTOF2024",
                "available": available_models
            },
            "analytics": {
                "enabled": False,
                "google-analytics": {"enabled": False, "tag": ""},
                "piwik-pro": {"enabled": False, "tag": ""}
            }
        }
        
        st.session_state.settings = default_settings
        st.success("🔄 Settings reset to defaults!")
        st.rerun()

with col3:
    st.info("💡 **Note**: Some settings require an application restart to take full effect.")

# Current Configuration Preview
st.markdown("---")  # Separator
st.header("📋 Current Configuration")

with st.expander("View Current Settings JSON"):
    st.json(st.session_state.settings)

# Resource Usage (if available)
try:
    import psutil
    st.markdown("---")  # Separator
    st.header("📊 System Resources")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cpu_percent = psutil.cpu_percent(interval=1)
        st.metric("CPU Usage", f"{cpu_percent:.1f}%")
    
    with col2:
        memory = psutil.virtual_memory()
        st.metric("Memory Usage", f"{memory.percent:.1f}%")
    
    with col3:
        disk = psutil.disk_usage('/')
        st.metric("Disk Usage", f"{disk.percent:.1f}%")
    
    with col4:
        st.metric("Available Memory", f"{memory.available / 1024**3:.1f} GB")

except ImportError:
    st.info("💡 Install `psutil` for system resource monitoring: `pip install psutil`")