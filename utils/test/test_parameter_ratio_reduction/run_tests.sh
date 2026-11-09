#!/bin/sh
# Production generic/specification/adapter tests; not a full SUMMA/MPI/TOML test.
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/../../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
cd "$work"
fc=${FC:-gfortran}
flags="${FFLAGS:--O0 -g -fcheck=all -fbacktrace} -ffree-line-length-none -cpp"
src="$root/build/source"
echo "Compiler: $($fc --version | head -1)"
echo "Flags: $flags"
# Only globalData and summa_type holders are replaced. Routines are production code.
# shellcheck disable=SC2086
$fc $flags -c \
  "$src/numrec/nr_type.f90" \
  "$src/dshare/var_lookup.f90" \
  "$src/calibration/calibration_types.f90" \
  "$src/dshare/data_types.f90" \
  "$here/test_context.f90" \
  "$src/calibration/parameter_search.f90" \
  "$src/calibration/parameter_fields.f90" \
  "$src/calibration/parameter_constraints.f90" \
  "$src/calibration/model_calibration.f90" \
  "$src/calibration/summa_parameter_fields.f90" \
  "$src/dshare/get_ixname.f90" \
  "$src/engine/apply_parameter_overrides.f90" \
  "$src/calibration/summa_parameter_spec.f90" \
  "$here/test_ratio_reduction.f90"
# shellcheck disable=SC2086
$fc $flags ./*.o -o test_ratios
./test_ratios
