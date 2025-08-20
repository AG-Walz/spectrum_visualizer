#!/bin/bash

# immunuXplore Application Launcher
# Simplified run script for the immunoXplore Streamlit application

echo "Starting immunuXplore..."

# Check if we're in the right directory
if [ ! -f "app.py" ]; then
    echo "Error: app.py not found. Please run this script from the immunoexplore directory."
    exit 1
fi

# Set environment variables
export BOKEH_ALLOW_WS_ORIGIN="*"

# Default port
PORT=${1:-8501}

echo "Launching immunuXplore on port $PORT..."
echo "Open your browser to: http://localhost:$PORT"

# Run the application
streamlit run app.py --server.port=$PORT --server.address=0.0.0.0