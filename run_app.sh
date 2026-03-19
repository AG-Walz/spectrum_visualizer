#!/bin/bash

# MSpecViz Application Launcher

echo "Starting MSpecViz..."

# Check if in the right directory
if [ ! -f "app.py" ]; then
    echo "Error: app.py not found. Please run this script from the MSpecViz directory."
    exit 1
fi

# Default port
PORT=${1:-8501}

echo "Launching MSpecViz on port $PORT..."
echo "Open your browser to: http://localhost:$PORT"

# Run the application
conda run --no-capture-output -n mspecviz \
    streamlit run app.py --server.port=$PORT

# NOTE: After pressing Ctrl+C you may see SQLite/SQLAlchemy tracebacks like:
#   "SQLite objects created in a thread can only be used in that same thread"
#   "Exception during reset or similar"
#
# These are EXPECTED and HARMLESS. They come from Streamlit's internal session
# state storage (SQLAlchemy + SQLite). On shutdown, Python's garbage collector
# runs in a different thread than the one that created the SQLite connections,
# which SQLite forbids — causing this noisy but inconsequential teardown error.
# No data is lost. You can safely ignore these messages.
