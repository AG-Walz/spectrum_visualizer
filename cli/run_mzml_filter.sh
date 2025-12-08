#!/bin/bash

# Launch mzML Filter Workspace
# This script starts the Streamlit interface for filtering mzML files based on peptide matches

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKSPACE_FILE="$SCRIPT_DIR/../workspaces/dataset_filter.py"

echo "🔬 Starting Dataset Filter Workspace..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

streamlit run "$WORKSPACE_FILE" \
    --server.port 8502 \
    --server.headless true \
    --browser.gatherUsageStats false

