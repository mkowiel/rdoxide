#!/bin/sh
# Regenerates docs/bench/upstream_speed.html with all sections; see tools/regen_speed_report.py.
#   tools/regen_speed_report.sh /path/to/rdkit-venv/bin/python
set -e
cd "$(dirname "$0")/.."
[ -n "$1" ] || { echo "usage: $0 RDKIT_PYTHON [extra regen_speed_report.py options]" >&2; exit 2; }
rd=$1; shift
exec python tools/regen_speed_report.py --rdkit-python "$rd" "$@"
