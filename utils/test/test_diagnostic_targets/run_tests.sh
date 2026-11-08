#!/usr/bin/env bash
# Focused tests of production metric/scalarizer/checker routines. Not a SUMMA build.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
BUILD=$(mktemp -d "${TMPDIR:-/tmp}/summa-diagnostics.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
FC=${FC:-gfortran}
# IEEE checks require normal floating-point semantics; do not use -ffast-math.
read -r -a FLAGS <<< "${FFLAGS:--O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -ffree-line-length-none}"
python3 "$HERE/extract_routines.py" "$ROOT" "$BUILD"
cd "$BUILD"
"$FC" "${FLAGS[@]}" -c "$ROOT/build/source/numrec/nr_type.f90" \
  "$ROOT/build/source/objfunc/metrics.f90" \
  "$ROOT/build/source/calibration/calibration_types.f90" \
  scalarizer_under_test.f90 settings_under_test.f90
"$FC" "${FLAGS[@]}" "$HERE/test_core.f90" ./*.o -o test_core
./test_core
if [[ -f "$HERE/scoring_stubs.f90" ]]; then
  "$FC" "${FLAGS[@]}" -c "$HERE/scoring_stubs.f90" scorer_under_test.f90
  "$FC" "${FLAGS[@]}" "$HERE/test_scoring.f90" ./*.o -o test_scoring
  ./test_scoring
fi
