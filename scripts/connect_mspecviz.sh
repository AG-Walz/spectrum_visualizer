#!/usr/bin/env bash
# Connect to MSpecViz via SSH tunnel through a jumphost.
#
# Usage (from your LOCAL machine):
#   ./connect_mspecviz.sh                        # defaults
#   ./connect_mspecviz.sh myuser@jumphost myuser@server
#   LOCAL_PORT=9501 ./connect_mspecviz.sh        # use a different local port
#
# Prerequisites:
#   - SSH access to the jumphost and the target server
#   - SSH config entries (recommended) or full user@host strings

set -euo pipefail

JUMPHOST="${1:-jumphost}"
SERVER="${2:-peptimine_highmem1}"
LOCAL_PORT="${LOCAL_PORT:-8501}"
REMOTE_PORT=8501

echo "Forwarding localhost:${LOCAL_PORT} -> ${SERVER}:${REMOTE_PORT} via ${JUMPHOST}"
echo "Open http://localhost:${LOCAL_PORT} in your browser"
echo "Press Ctrl+C to disconnect"
echo ""

# Open browser after a short delay (works on macOS and Linux)
(
  sleep 2
  if command -v xdg-open &>/dev/null; then
    xdg-open "http://localhost:${LOCAL_PORT}" 2>/dev/null &
  elif command -v open &>/dev/null; then
    open "http://localhost:${LOCAL_PORT}" 2>/dev/null &
  fi
) &

# Establish the tunnel (-N = no remote command, just forward)
ssh -N -L "${LOCAL_PORT}:localhost:${REMOTE_PORT}" -J "${JUMPHOST}" "${SERVER}"
