#!/bin/sh
# builds coordgenlibs 3.0.2 (without the Maestro reader) as objects in OBJDIR for the oracle drivers:
#   tools/coordgen_drv/build.sh SRCDIR OBJDIR
set -e
src=$1; obj=$2
mkdir -p "$obj"
for f in "$src"/*.cpp; do
  b=$(basename "$f" .cpp)
  c++ -std=c++17 -O1 -ffp-contract=off -w -c "$f" -o "$obj/$b.o" &
done
wait
