#!/bin/sh
# Optional parity check against the PREVIOUS BLOCK PROTOTYPE, not live SUMMA.
set -eu
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/extracted/summa_feasibility_design" >&2
  exit 2
fi
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ref=$(CDPATH= cd -- "$1" && pwd)/build/source
if [ -d "$here/support" ]; then
  nr="$here/support/nr_type.f90"
  search="$here/support/parameter_search.f90"
else
  root=$(CDPATH= cd -- "$here/../../.." && pwd)
  nr="$root/build/source/numrec/nr_type.f90"
  search="$root/build/source/calibration/parameter_search.f90"
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
fc=${FC:-gfortran}
flags="-cpp -O0 -g -fcheck=all -ffree-line-length-none"
mkdir "$work/old" "$work/new"
cd "$work/old"
# shellcheck disable=SC2086
$fc $flags "$ref/numrec/nr_type.f90" "$ref/calibration/parameter_search.f90" \
  "$ref/calibration/parameter_fields.f90" "$ref/calibration/parameter_constraints.f90" \
  "$ref/calibration/model_calibration.f90" "$here/tests/parity_fixture.F90" -o parity
./parity > "$work/old.txt"
cd "$work/new"
# shellcheck disable=SC2086
$fc $flags -DNESTED "$nr" "$search" \
  "$here/src/parameter_fields.f90" "$here/src/parameter_constraints.f90" \
  "$here/src/model_calibration.f90" "$here/tests/parity_fixture.F90" -o parity
./parity > "$work/new.txt"
# Compare every decision exactly, and field sums with a relative rounding allowance
# because nested versus flattened reductions can change floating addition order.
python3 - "$work/old.txt" "$work/new.txt" <<'PY'
import math
import sys
from pathlib import Path
old = [s.split() for s in Path(sys.argv[1]).read_text().splitlines()]
new = [s.split() for s in Path(sys.argv[2]).read_text().splitlines()]
assert len(old) == len(new)
for i, (a, b) in enumerate(zip(old, new), 1):
    assert a[0] == b[0], ('feasibility mismatch', i, a, b)
    assert math.isclose(float(a[1]), float(b[1]), rel_tol=1e-14, abs_tol=1e-14), ('field mismatch', i, a, b)
print(f'PASS: block/nested parity ({len(old)} candidate decisions and field sums).')
PY
