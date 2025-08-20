import shutil
from pathlib import Path

import streamlit as st
import pandas as pd

import sys

# Add src to path for imports
src_dir = Path(__file__).parent.parent / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from src.common.common import reset_directory