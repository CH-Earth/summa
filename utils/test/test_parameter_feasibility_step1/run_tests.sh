#!/bin/sh
# Regression tests for the PRODUCTION generic feasibility and SUMMA adapter modules.
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/../../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
cd "$work"
fc=${FC:-gfortran}
flags="${FFLAGS:--O0 -g -fcheck=all -fbacktrace -Wall -Wextra} -ffree-line-length-none"
echo "Compiler: $($fc --version | head -1)"
echo "Flags: $flags"
# Build the real production modules, not copied versions under utils/.
# shellcheck disable=SC2086
$fc $flags -c "$root/build/source/numrec/nr_type.f90" \
  "$root/build/source/calibration/parameter_search.f90" \
  "$root/build/source/calibration/parameter_fields.f90" \
  "$root/build/source/calibration/parameter_constraints.f90" \
  "$root/build/source/calibration/model_calibration.f90" \
  "$root/build/source/calibration/nsga2.f90" \
  "$here/tests/test_support.f90"
common="nr_type.o parameter_search.o parameter_fields.o parameter_constraints.o model_calibration.o nsga2.o test_support.o"
# shellcheck disable=SC2086
$fc $flags $common "$here/tests/test_core.f90" -o test_core
./test_core
# shellcheck disable=SC2086
$fc $flags $common "$here/tests/test_optimizers.f90" -o test_optimizers
./test_optimizers
# Real SUMMA types, with the unrelated globalData module replaced by a small stub.
# shellcheck disable=SC2086
$fc $flags -cpp -c "$root/build/source/dshare/var_lookup.f90" \
  "$root/build/source/calibration/calibration_types.f90" \
  "$root/build/source/dshare/data_types.f90" \
  "$here/tests/global_data_stub.f90" \
  "$root/build/source/calibration/summa_parameter_fields.f90"
# shellcheck disable=SC2086
$fc $flags $common var_lookup.o calibration_types.o data_types.o global_data_stub.o \
  summa_parameter_fields.o "$here/tests/test_adapter.f90" -o test_adapter
./test_adapter
