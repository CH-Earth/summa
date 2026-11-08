#!/usr/bin/env bash
# Tests the production metric and extracted production scorer in a light harness.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
DIAG="$ROOT/utils/test/test_diagnostic_targets"
BUILD=$(mktemp -d "${TMPDIR:-/tmp}/summa-boxcox.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
FC=${FC:-gfortran}
# Do not use -ffast-math: these tests depend on IEEE finite/NaN semantics.
read -r -a FLAGS <<< "${FFLAGS:--O0 -g -fcheck=all,no-array-temps -fbacktrace -ffpe-trap=invalid,zero,overflow -ffree-line-length-none}"
python3 "$HERE/check_interfaces.py" "$ROOT"
python3 "$DIAG/extract_routines.py" "$ROOT" "$BUILD"
cd "$BUILD"
"$FC" "${FLAGS[@]}" -c "$ROOT/build/source/numrec/nr_type.f90" \
  "$ROOT/build/source/objfunc/metrics.f90" \
  "$ROOT/build/source/calibration/calibration_types.f90" scalarizer_under_test.f90
"$FC" "${FLAGS[@]}" "$HERE/test_metrics.f90" ./*.o -o test_metrics
./test_metrics
"$FC" "${FLAGS[@]}" -c "$DIAG/scoring_stubs.f90" scorer_under_test.f90
"$FC" "${FLAGS[@]}" "$HERE/test_scoring.f90" ./*.o -o test_scoring
./test_scoring
python3 "$HERE/extract_config.py" "$ROOT" "$BUILD"
"$FC" "${FLAGS[@]}" -c "$HERE/toml_stubs.f90" parsers_under_test.f90 legacy_under_test.f90
"$FC" "${FLAGS[@]}" "$HERE/test_config.f90" ./*.o -o test_config
./test_config
# Isolate the real evaluation-writer module from the scorer harness substitute.
mkdir "$BUILD/metadata"
cd "$BUILD/metadata"
python3 "$HERE/extract_output.py" "$ROOT" "$BUILD/metadata"
"$FC" "${FLAGS[@]}" -c "$ROOT/build/source/numrec/nr_type.f90" \
  "$ROOT/build/source/objfunc/metrics.f90" \
  "$ROOT/build/source/calibration/calibration_types.f90" \
  "$ROOT/build/source/calibration/parameter_search.f90" "$HERE/netcdf_stubs.f90" \
  "$ROOT/build/source/objfunc/write_evaluation.f90" calibration_output_under_test.f90
"$FC" "${FLAGS[@]}" "$HERE/test_output.f90" ./*.o -o test_output
./test_output
