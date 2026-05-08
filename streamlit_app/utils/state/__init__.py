"""Session state, workspaces, and page setup for the Streamlit app.

Re-exports every public name from the four submodules so callers keep using
``from streamlit_app.utils.state import ...`` without caring which file a
name lives in.
"""
from streamlit_app.utils.state.state_core import (
    OS_PLATFORM,
    PROJECT_ROOT,
    STATE_PREFIX,
    _persisted_key,
    get_param,
    get_state_key,
)
from streamlit_app.utils.state.state_workspace import (
    is_safe_workspace_name,
)
from streamlit_app.utils.state.state_params import (
    load_params,
    save_params,
)
from streamlit_app.utils.state.state_sidebar import (
    page_setup,
    render_sidebar,
)

__all__ = [
    # core
    "OS_PLATFORM",
    "PROJECT_ROOT",
    "STATE_PREFIX",
    "_persisted_key",
    "get_param",
    "get_state_key",
    # workspace
    "is_safe_workspace_name",
    # params / setup
    "load_params",
    "page_setup",
    "render_sidebar",
    "save_params",
]
