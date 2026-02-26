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
    streamlit run app.py --server.port=$PORT --server.address=0.0.0.0
