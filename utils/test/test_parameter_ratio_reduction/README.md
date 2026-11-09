# Explicit local reference reduction for chained parameter ratios

## Scope

This change was prepared against the user's uploaded `summa_calibration_source.tar`.
It is not a patch against a live checkout; use `git apply --check` first.
No MPI scheduler, model storage layout, optimizer or model-physics code changes.

## Semantics

`reference_reduction` is an optional entry in an inline parameter specification
under `[calibration.parameter_ratios]`. It defaults to `"none"`.

- `none`: retain existing elementwise construction and equal local length checks.
- `min`: the dependent must be a local scalar; at each horizontal location compute
  `dependent(h,1) = ratio * minval(reference_effective(h)%dat(:))`.

The minimum is over local components/layers, NEVER over horizontal locations.
It is taken BEFORE multiplication by the ratio. The reference is effective, so
upstream calibrated multipliers/ratios have already been applied. Each candidate
starts from immutable initialized fields; the downstream scalar is not expanded.
SUMMA additionally requires the model-parameter target to have native `scalarv`
metadata (a one-layer `parSoil` parameter is not treated as a native scalar).

For positive saturation and ratios strictly between zero and one, the example
ensures wilting < transpiration < saturation in every layer. It does not by itself
ensure residual moisture < wilting or a configured physical minimum gap. Retain
ordered constraints for those conditions. The patch removes the specification
builder's obsolete bans on chained ratios and ratio targets in ordered chains;
generic dependency sorting, cycle detection and physical constraint checks remain.

Constraint endpoint ranges use the same reduced reference: the extrema across
LOCAL minima, not extrema of all source layer values. The calibration NetCDF
metadata records `reference_reduction` next to `reference_name`.

This is an opt-in modeling choice, not an assumption that a minimum is the best
scientific parameterization for every model or application. No mean, root-zone
weighting, layer selection or scalar-to-layer ratio broadcasting is introduced.

## Run

From the repository root:

```sh
./utils/test/test_parameter_ratio_reduction/run_tests.sh
FFLAGS='-O2 -g -fcheck=all -fbacktrace' \
  ./utils/test/test_parameter_ratio_reduction/run_tests.sh
```

The tests compile the real generic field/constraint modules, the SUMMA
specification builder, native data types, name lookup, and field adapter.
Only `globalData` and `summa_type` holders are replaced by minimal test modules.
The test context file must NEVER be compiled into the production executable.

Coverage includes reverse registry/dependency order, heterogeneous layer counts,
local minima, repeated samples with no accumulated scaling, native scalar
write-back, unchanged fixed fields, fixed references, signed reductions,
reduction-aware gap extrema, cycles, self-reference, invalid reduction metadata,
and incompatible shapes/native metadata.

## Verification in this workspace

GNU Fortran 14.2.0, debug and optimized builds with `-fcheck=all`:

- Existing nested feasibility core: 4,190 assertions passed.
- Existing random/DDS/NSGA-II integration: 431 assertions passed.
- Existing SUMMA native adapter: 68 assertions passed.
- New local-minimum/chain/specification/adapter suite: 118 assertions passed.

The existing tests were restored from the earlier review bundle because the
uploaded archive contains production source but not `utils/test`.

The example TOML fragment was syntax-checked with Python `tomllib` only.
The new Fortran TOML parsing branch and NetCDF metadata-writing branch were
reviewed but not compiled/executed here: this source-only archive does not include
the required external dependencies. A full SUMMA build, actual TOML parsing,
NetCDF readback and short Athabasca MPI run remain required before committing.

The change covers the current field-based calibration path. The obsolete
`build_summa_parameter_overrides` helper is not used by that path and has not been
extended into a layered-field or chained-ratio API.
