#!/bin/sh
# Unit tests of the calibration modules, built from the production sources (no SUMMA run, MPI, TOML or NetCDF).
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
src=$(CDPATH= cd -- "$here/../../../build/source" && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/summa-calib-units.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cd "$work"
fc=${FC:-gfortran}
# No -ffast-math: the metric tests rely on IEEE NaN and overflow behaviour.
flags="${FFLAGS:--O0 -g -fcheck=all -fbacktrace} -ffree-line-length-none -cpp"
echo "Compiler: $($fc --version | head -1)"
echo "Flags: $flags"

# Production modules; test_context.f90 stands in for globalData and summa_type only.
# shellcheck disable=SC2086
$fc $flags -c \
  "$src/numrec/nr_type.f90" \
  "$src/dshare/var_lookup.f90" \
  "$src/calibration/calibration_types.f90" \
  "$src/dshare/data_types.f90" \
  "$here/test_context.f90" \
  "$src/objfunc/metrics.f90" \
  "$src/calibration/parameter_search.f90" \
  "$src/calibration/parameter_fields.f90" \
  "$src/calibration/parameter_constraints.f90" \
  "$src/calibration/model_calibration.f90" \
  "$src/calibration/nsga2.f90" \
  "$src/calibration/summa_parameter_fields.f90" \
  "$src/dshare/get_ixname.f90" \
  "$src/engine/apply_parameter_overrides.f90" \
  "$src/calibration/summa_parameter_spec.f90" \
  "$here/test_support.f90"
objs=$(ls ./*.o)

for t in test_metrics test_fields test_optimizers test_adapter test_ratio_reduction; do
  # shellcheck disable=SC2086
  $fc $flags $objs "$here/$t.f90" -o "$t"
  "./$t"
done
