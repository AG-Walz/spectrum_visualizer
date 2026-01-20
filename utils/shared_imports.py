"""
Centralized package imports - Loads all packages once at startup.
Imports packages to prevent redundant loading when navigating between pages.
"""

# Core scientific packages
import pandas as pd
import numpy as np

# MS/MS analysis packages (heavy imports)
import pyopenms as oms
from psm_utils.io import read_file
from psm_utils import Peptidoform, PSM, PSMList

# Standard library imports
import argparse
import warnings
import tempfile
import os
import sys
import logging
import json
import gc
import shutil
import uuid
import time
import zipfile
import socket
import datetime
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any, Union, Set, Iterable
from collections import defaultdict
import itertools
import re
from io import StringIO, BytesIO
import io

# Optional packages
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    psutil = None

# Bokeh components (import bokeh first, then streamlit-bokeh)
try:
    from bokeh.plotting import figure
    from bokeh.embed import file_html
    from bokeh.resources import CDN
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False
    figure = None
    file_html = None
    CDN = None

# Streamlit-Bokeh integration (optional)
try:
    from streamlit_bokeh import st_bokeh
    STREAMLIT_BOKEH_AVAILABLE = True
except ImportError:
    STREAMLIT_BOKEH_AVAILABLE = False
    st_bokeh = None

# Tabulate
try:
    from tabulate import tabulate
    TABULATE_AVAILABLE = True
except ImportError:
    TABULATE_AVAILABLE = False
    tabulate = None

# Handle NumPy compatibility
warnings.filterwarnings("ignore", message=".*np.bool8.*", category=DeprecationWarning)
if not hasattr(np, "bool8"):
    np.bool8 = np.bool_

# Configure logging (do this once)
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# Configure environment
os.environ.setdefault("BOKEH_ALLOW_WS_ORIGIN", "*")

# Add core and utils directories to path
utils_dir = Path(__file__).parent
core_dir = Path(__file__).parent.parent / "core"
if str(utils_dir) not in sys.path:
    sys.path.insert(0, str(utils_dir))
if str(core_dir) not in sys.path:
    sys.path.insert(0, str(core_dir))

# Import common utilities from utils/helpers.py
try:
    from helpers import (
        create_file_options_with_sizes, map_display_names_to_paths,
        extract_filename, extract_filenames, is_empty_collection,
        safe_concat_dataframes, show_status_message, validate_file_list,
        create_validation_info_dict
    )
except ImportError:
    # Fallback for when imported from different context
    create_file_options_with_sizes = None
    map_display_names_to_paths = None
    extract_filename = None
    extract_filenames = None
    is_empty_collection = None
    safe_concat_dataframes = None
    show_status_message = None
    validate_file_list = None
    create_validation_info_dict = None

# Spectrum tools (Moved to individual files to avoid circular imports)
SPECTRUM_TOOLS_AVAILABLE = False

# Export all commonly used items
__all__ = [
    # Packages
    'pd', 'np', 'oms', 'read_file', 'Peptidoform', 'PSM', 'PSMList',
    'argparse', 'warnings', 'tempfile', 'os', 'sys', 'logging', 'json', 'gc',
    'shutil', 'uuid', 'time', 'zipfile', 'socket', 'datetime',
    'Path', 'Optional', 'Tuple', 'List', 'Dict', 'Any', 'Union', 'Set', 'Iterable',
    'defaultdict', 'itertools', 're', 'io', 'StringIO', 'BytesIO',
    # Optional packages
    'psutil', 'PSUTIL_AVAILABLE',
    'file_html', 'CDN', 'figure', 'BOKEH_AVAILABLE',
    'st_bokeh', 'STREAMLIT_BOKEH_AVAILABLE',
    'tabulate', 'TABULATE_AVAILABLE',
    # Logger
    'logger',
    # Common utilities
    'create_file_options_with_sizes', 'map_display_names_to_paths',
    'extract_filename', 'extract_filenames', 'is_empty_collection',
    'safe_concat_dataframes', 'show_status_message', 'validate_file_list',
    'create_validation_info_dict',
    # Spectrum tools
    'SPECTRUM_TOOLS_AVAILABLE',
]