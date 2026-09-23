#!/bin/sh
# Start the EASLEE web GUI: http://localhost:8080
cd "$(dirname "$0")" || exit 1
exec julia --threads=auto --project=webgui webgui/easlee_app.jl "$@"
