#!/usr/bin/env bash
# Re-sync the vendored solver and analysis modules from pgcp-experiments
# whenever the upstream files change. The webapp keeps its own copies under
# webapp/server/ so the container can be built without the pgcp-experiments
# tree alongside; this script is the one place that knows where the upstream
# files live.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../../pgcp-experiments/scripts"
DST="$HERE/../server"
for f in bdd_solver_triangle.jl bdd_analysis.jl; do
    if ! [ -f "$SRC/$f" ]; then
        echo "missing $SRC/$f — is the pgcp-experiments tree present?" >&2
        exit 1
    fi
    cp -v "$SRC/$f" "$DST/$f"
done
echo "synced into $DST"
