#!/usr/bin/env bash
# Build Cyndi (BMC Bioinformatics 2009) with MSYS2/MinGW g++ on Windows.
# Usage: bash build.sh [output-name]
# Notes:
#  - The MSYS2 g++ at C:/Users/Administrator/tools/mingw64 does NOT auto-include
#    its own sysroot headers; -idirafter must point at <mingw>/include explicitly.
#  - Original code has one portability fix: MMFF94StretchBend.cpp uses
#    `dynamic_cast<MMFF94*>(...) == 0` instead of `== false` (GCC 16 rejects the latter).
#  - Cyndi.cpp OptimizeConformer branch: restored FF setup()/cgm.setup() calls that were
#    commented out in the original (they are required for CG minimization to actually run).
set -e
cd "$(dirname "$0")"
export PATH="/c/Users/Administrator/tools/mingw64/bin:$PATH"
OUT="${1:-build/cyndi.exe}"
MINGW='C:/Users/Administrator/tools/mingw64/include'

mkdir -p build
g++ -O3 -DWIN32 -DNDEBUG -D_CONSOLE \
    -I include \
    -idirafter "$MINGW" \
    src/*.cpp \
    -o "$OUT"
# The MSYS2 runtime DLLs are not on PATH by default; copy them next to the
# exe so it runs from any shell (Windows loads DLLs from the exe's dir first).
MINGWBIN="/c/Users/Administrator/tools/mingw64/bin"
for dll in libstdc++-6.dll libgcc_s_seh-1.dll libwinpthread-1.dll; do
    cp -f "$MINGWBIN/$dll" "$(dirname "$OUT")/" 2>/dev/null || true
done
echo "Built: $OUT"
echo "Note: define TORSION_LINEAR to disable circular torsion encoding (-DTORSION_LINEAR)"
