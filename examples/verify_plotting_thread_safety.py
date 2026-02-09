import sys
from pathlib import Path
import time
import pandas as pd
import numpy as np
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

# Setup paths
cwd = Path.cwd()
specviz_dir = Path("/mnt/volume/nacho/specviz/spectrum_visualizer")
sys.path.insert(0, str(specviz_dir / "utils"))
sys.path.insert(0, str(specviz_dir / "core"))

# Mock streamlit
from unittest.mock import MagicMock
if 'streamlit' not in sys.modules:
    sys.modules['streamlit'] = MagicMock()
def noop_decorator(func=None, **kwargs):
    if func is None: return lambda f: f
    return func
sys.modules['streamlit'].cache_data = noop_decorator
sys.modules['streamlit'].cache_resource = noop_decorator

try:
    from spectrum_engine import create_matplotlib_mirror_plot
    print("✅ Imported create_matplotlib_mirror_plot")
except ImportError as e:
    print(f"❌ Failed to import spectrum_engine: {e}")
    sys.exit(1)

# Generate dummy data
def generate_dummy_spectrum(n_peaks=100):
    mz = np.linspace(100, 1000, n_peaks)
    intensity = np.random.rand(n_peaks)
    annotation = [""] * n_peaks
    # Annotate some peaks
    for i in range(0, n_peaks, 10):
        annotation[i] = f"y{i//10}"
    return pd.DataFrame({"mz": mz, "intensity": intensity, "annotation": annotation})

print("Generating dummy data...")
upper_df = generate_dummy_spectrum()
lower_df = generate_dummy_spectrum()

def test_single_plot(i):
    try:
        fig = create_matplotlib_mirror_plot(
            upper_df, lower_df, 
            title=f"Test Plot {i}",
            r2_matching=0.95, r2_binned=0.92
        )
        return True
    except Exception as e:
        print(f"❌ Error in plot {i}: {e}")
        return False

# Sequential Test
n_plots = 20
print(f"\nRunning sequential test ({n_plots} plots)...")
start = time.time()
for i in range(n_plots):
    test_single_plot(i)
seq_time = time.time() - start
print(f"Sequential time: {seq_time:.4f}s ({seq_time/n_plots:.4f}s/plot)")

# Parallel Test
print(f"\nRunning parallel test ({n_plots} plots)...")
start = time.time()
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(test_single_plot, i) for i in range(n_plots)]
    results = [f.result() for f in concurrent.futures.as_completed(futures)]
par_time = time.time() - start
print(f"Parallel time: {par_time:.4f}s ({par_time/n_plots:.4f}s/plot)")

speedup = seq_time / par_time
print(f"\nSpeedup: {speedup:.2f}x")

if all(results):
    print("✅ Parallel execution successful (thread-safe)")
else:
    print("❌ Parallel execution failed")
