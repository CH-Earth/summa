# test_calibration_units

Unit tests of the calibration modules. They compile the production sources in `build/source/`
directly, with gfortran alone: no SUMMA build, MPI, TOML or NetCDF. About 15 s.

```sh
./utils/test/test_calibration_units/run_tests.sh
FFLAGS='-O2 -g -fcheck=all -fbacktrace' ./utils/test/test_calibration_units/run_tests.sh
```

The build happens in a temporary directory, so the checkout gets no `.o` or `.mod` files.
`FC` picks the compiler. Do not add `-ffast-math`: the metric tests depend on IEEE NaN and
overflow behaviour.

| Test | Covers |
|---|---|
| `test_metrics.f90` | `metrics.f90`: every metric and transformation, the Box-Cox exponent (default, formula, limits near zero, overflow guards) and the undefined results that come back as NaN with `metric_undefined` |
| `test_fields.f90` | `parameter_fields`, `parameter_constraints`, `model_calibration`: nested HRU/layer fields, spatial multipliers, ratios, ordered chains and the feasibility check |
| `test_optimizers.f90` | DDS, random sampling and NSGA-II with the same feasibility check, redrawing infeasible candidates |
| `test_adapter.f90` | `summa_parameter_fields`: reading SUMMA parameter structures into fields and writing them back |
| `test_ratio_reduction.f90` | `summa_parameter_spec` and the adapter together: chained ratios with `reference_reduction = "min"`, cycles, bad metadata, shapes |

`test_context.f90` stands in for `globalData` and `summa_type`, and `test_support.f90` holds shared
check helpers. Never add either to a SUMMA build.

What these do not cover: TOML parsing, scoring a simulation against observations, and the calibration
NetCDF output. `../test_calibration/test_calibration_bow_pareto.sh` covers those end to end,
including a zero-weight Box-Cox diagnostic target on the DDS run.
