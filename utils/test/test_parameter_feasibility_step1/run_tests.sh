#!/bin/sh
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
# Bundle: use frozen, unmodified uploaded dependencies. Additive-checkout patch:
# use the checkout's existing shared types/search module without editing them.
if [ -d "$here/support" ]; then
  nr="$here/support/nr_type.f90"
  search="$here/support/parameter_search.f90"
  lookup="$here/support/var_lookup.f90"
  types="$here/support/data_types.f90"
  ctypes="$here/support/calibration_types.f90"
else
  root=$(CDPATH= cd -- "$here/../../.." && pwd)
  nr="$root/build/source/numrec/nr_type.f90"
  search="$root/build/source/calibration/parameter_search.f90"
  lookup="$root/build/source/dshare/var_lookup.f90"
  types="$root/build/source/dshare/data_types.f90"
  ctypes="$root/build/source/calibration/calibration_types.f90"
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
cd "$work"
fc=${FC:-gfortran}
# IEEE overflow/NaN candidates are checked and rejected, not trapped. Do not use
# -ffast-math or -ffinite-math-only: they invalidate IEEE finite-value checking.
flags="${FFLAGS:--O0 -g -fcheck=all -fbacktrace -Wall -Wextra} -ffree-line-length-none"
echo "Compiler: $($fc --version | head -1)"
echo "Flags: $flags"
# Intentional word splitting allows normal FC/FFLAGS overrides.
# shellcheck disable=SC2086
$fc $flags -c "$nr" "$search" \
  "$here/src/parameter_fields.f90" "$here/src/parameter_constraints.f90" \
  "$here/src/model_calibration.f90" "$here/tests/test_support.f90"
# shellcheck disable=SC2086
$fc $flags nr_type.o parameter_search.o parameter_fields.o parameter_constraints.o \
  model_calibration.o test_support.o "$here/tests/test_core.f90" -o test_core
./test_core
if [ -f "$here/tests/test_adapter.f90" ]; then
  # Uses real SUMMA types; only globalData is substituted with its small test holder.
  # shellcheck disable=SC2086
  $fc $flags -cpp -c "$lookup" "$ctypes" \
    "$types" "$here/tests/global_data_stub.f90" \
    "$here/src/summa_parameter_fields.f90"
  # shellcheck disable=SC2086
  $fc $flags nr_type.o parameter_search.o parameter_fields.o parameter_constraints.o \
    model_calibration.o test_support.o var_lookup.o calibration_types.o data_types.o \
    global_data_stub.o summa_parameter_fields.o "$here/tests/test_adapter.f90" -o test_adapter
  ./test_adapter
fi
