# Zero-weight diagnostic targets for SUMMA DDS

## Behavior

All targets are evaluated from the same model run and returned through the existing
MPI objective vector. DDS uses only positive-weight targets. Zero-weight targets
are evaluated and saved but do not enter elite selection.

- `scalarize_objectives` skips a zero weight BEFORE accessing its metric.
- `check_search_settings` requires finite nonnegative weights and at least one
  positive-weight target for optimization.
- NSGA-II currently rejects any zero-weight target, rather than silently letting
  a diagnostic affect dominance or crowding. Positive-weight NSGA-II settings
  retain their existing interpretation; magnitudes of weights are not used there.
- `target_weight(target)` is stored alongside `target_name` in calibration NetCDF.
- A mathematically undefined diagnostic is NaN in memory, logged, and stored as
  the objective variable's explicit NetCDF `_FillValue`. `trial_failed` remains 0.
- An undefined ACTIVE metric marks that trial as failed. This replaces the old
  finite +/-1e6 placeholders for undefined efficiencies/errors. Valid metric
  formulas are unchanged; genuine very poor but finite scores are not clipped.
- Configuration, missing-file, variable, alignment, and output errors are NOT
  silently ignored for diagnostic targets. Model-physics failure handling is unchanged.

The scope is diagnostic-only scoring in DDS, not a new low-flow metric suite or
active-objective subset support for NSGA-II. NSE-log uses the existing offset
`mean(observed valid pairs)/100`. Invalid transform domains and constant-series
undefined efficiencies are detected explicitly.

## Source basis

The production patch was prepared against `summa_calibration_source.tar`, uploaded
November 8, 2026, with `diff(20261108-100534).txt` applied first. It changes five
production files, and includes this test directory. It does not change parameter
constraints, physics, MPI messages, or the TOML target schema. Your local checkout
may differ: use `git apply --check` before applying it. No actual project commit
hash is claimed for the provided source archive.

## Apply and build

From the repository root:

```sh
git apply --check /path/to/summa_zero_weight_diagnostics.patch
git apply /path/to/summa_zero_weight_diagnostics.patch
git diff --check
bash utils/test/test_diagnostic_targets/run_tests.sh
```

Then rebuild the normal SUMMA executables using your existing build configuration.
Do not build these IEEE-validity checks with `-ffast-math` / finite-math-only
assumptions. Do not overwrite full source files from another revision to avoid
an application conflict; inspect the affected hunks instead.

## Configuration

Append the two blocks in `diagnostic_targets.toml` to the current exp50 template
(or edit an existing explicit target list). Preserve `algorithm = "dds"` and the
existing parameter sections. Include the ordinary streamflow-KGE target explicitly:
the parser does not append the legacy implicit target when explicit targets exist.

`write_timeseries = false` can remain unchanged: every target value is written to
the calibration output file without full model time-series output. Retain the
small evaluation budget for a smoke test, then explicitly set the production
budget to 5000 when ready; the patch does not modify any evaluation budget.

## Tests included

`run_tests.sh` uses gfortran and Python 3. It compiles the actual `nr_type`,
`calibration_types` and `metrics` modules, and extracts the production scalarizer,
search-setting validator and `score_target` into a focused test harness.

The target-scoring harness substitutes observation I/O, time alignment and
standalone evaluation output. Thus it tests the production scoring/error-handling
logic, but is NOT a full SUMMA, MPI, or NetCDF integration test. Baseline removal
and seasonal balance calls are intentionally not exercised; the substitutes fail
if those paths are invoked accidentally.

Expected result:

```
PASS: metrics, scalarizer and optimization settings (48 assertions).
PASS: production target scorer with test I/O/alignment substitutes (16 assertions).
```

These 64 assertions cover valid log NSE including zero-flow observations, undefined
metrics, unknown metric/transform names, empty/missing pairs, scalarization of
zero-weight NaN/Inf, weight validation, the NSGA-II guard, diagnostic warnings,
active-target failure, and propagation of input/alignment/output errors.

An additional validation against the original archived metrics module checked
1500 valid scores (5 metrics x 6 transformations x 50 series pairs), all bit-for-bit
identical under gfortran at -O0. This comparison used the original source archive
and is not re-run by the bundled harness. It does not establish end-to-end model
or asynchronous-search equivalence.

## Required full-model checks on the working installation

First compare the primary KGE and DDS scalar for IDENTICAL parameter sets with
and without the zero-weight NSE-log target. Then run a small DDS calibration with
the example target list and check its actual NetCDF output:

```sh
python3 utils/test/test_diagnostic_targets/check_output.py /path/to/calibration.nc
```

The optional checker requires numpy and netCDF4. It validates target weights,
completed/failed trials, and finite active objectives. Missing diagnostic values
are reported, not treated as trial failures. It does not recompute the metrics.

A full model/MPI run and the real NetCDF-Fortran writer were not executed in the
patch-building environment. The focused tests do not substitute for that check.

## Plotting and experiment comparisons

Read target names and weights rather than assuming a single objective column.
Select the best parameter sample using ordinary KGE, then report the diagnostic
at THAT sample. Do not select a separate maximum-NSE-log sample for the primary
comparison. Do not reject a successful trial solely because a zero-weight target
is missing, and report the number of defined diagnostic values for each comparison.

Adding target-scoring work can change worker completion order in asynchronous DDS,
so an identical RNG seed is not a guarantee of an identical optimization trajectory.
Use fixed parameter evaluations or a single evaluation worker for strict objective
parity tests.
