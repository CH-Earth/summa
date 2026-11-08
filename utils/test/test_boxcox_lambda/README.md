# Configurable Box-Cox exponent for SUMMA calibration targets

## Scope and prerequisite

This is an incremental change **after** `summa_zero_weight_diagnostics.patch`.
It is based on the November 8 source archive, the subsequent scalar-ratio diff,
and the diagnostic-target patch already applied in this conversation.

`boxcox_lambda` is a fixed TOML setting for evaluating a target, not a physical
parameter sampled by DDS or NSGA-II. All existing parameter bounds, optimizer
settings and zero-weight diagnostic behavior remain unchanged.

## Configuration

A target can override the exponent:

```toml
[[calibration.target]]
name          = "streamflow_nse_boxcox"
variable      = "streamflow"
obs_path      = "{home}/data/camels-spat/observations/obs-daily/"
obs_file      = "{case_name}_daily_flow_observations.nc"
vname_obs     = "q_obs"
metric        = "nse"
obs_transform = "boxcox"
boxcox_lambda = 0.2
weight        = 0.0
```

The value 0.2 is an example, not a newly recommended experimental bound. Append
this to the existing explicit target list, retaining the ordinary KGE target
with weight 1.0 and any existing diagnostic targets.

An optional `boxcox_lambda` key in the **existing** `[calibration]` section sets
the default inherited by targets that omit it. A target-local value takes
precedence. When neither location supplies it, the exponent remains **0.25**.
The legacy implicit single-target configuration inherits the same default.
Use real-valued TOML literals such as `0.0`, `0.2` or `-0.5`.

Do not place this key under `[calibration.parameters]`. It is not an additional
calibrated model parameter. It has no effect on targets using another transform.
Nonfinite values (NaN or infinity) and wrong TOML value types are configuration
errors, including for a zero-weight diagnostic.

## Mathematical definition and compatibility

The transform remains **unshifted**:

    T(Q; lambda) = (Q**lambda - 1) / lambda   for lambda != 0
                  log(Q)                   for lambda == 0

* lambda > 0: Q = 0 is allowed, with T(0) = -1/lambda.
* lambda <= 0: Q must be strictly positive.
* Negative Q is rejected for every exponent; zeros are not silently discarded.
* No flow offset or rescaling has been added.

In particular, `obs_transform="boxcox", boxcox_lambda=0.0` is **not** the existing
`obs_transform="log"`, which uses log(Q + epsilon), epsilon = mean(observed Q)/100.
That existing log behavior is unchanged.

The historical exponent 0.25 retains the exact old power expression. For other
exponents, a small-argument expansion avoids loss of precision near lambda=0.
Overflow guards return an undefined metric rather than evaluating an invalid
intermediate power. The zero-weight patch then handles undefined diagnostics
as missing without discarding an otherwise valid trial. An undefined active
objective retains the failed-trial behavior from that patch.

## Fortran interfaces: error outputs last

Optional inputs precede the error outputs. For example:

```fortran
subroutine compute_metric(obs,sim,metric,transfo,objective,boxcox_lambda,err,message)
```

Every affected call uses keyword error arguments, including when the optional
lambda is omitted:

```fortran
! Keep the historical default of 0.25.
call compute_metric(obs,sim,metric,transfo,objective, &
                    err=err,message=message)

! Override lambda explicitly.
call compute_metric(obs,sim,metric,transfo,objective, &
                    boxcox_lambda=0.2_rkind,err=err,message=message)
```

The same order is used by `prepare_series`, `apply_transformation`,
`write_evaluation`, and `define_evaluation_target`. The last helper retains its
existing dummy name `ierr`, so its call ends with `ierr=err,message=cmessage`.
Optional lambda is forwarded by keyword through the helper chain; an omitted
lambda stays absent until the default is selected.

This is an intentional source-interface change: old positional calls of the form
`...,objective,err,message)` must be updated. All direct calls in the available
source archive and both supplied test suites have been updated. External callers
and tests outside that archive also need the keyword form. Rebuild dependent
Fortran modules; do not mix old objects/module files with the new interface.

No transformation arithmetic, default value, TOML setting or output metadata
semantics changed relative to the original configurable-lambda patch.

## Output metadata

The calibration NetCDF file gains `target_boxcox_lambda(target)` in the same
order as `target_name`. Box-Cox targets contain their effective exponent; other
targets contain `NF90_FILL_DOUBLE`, declared in `_FillValue`. Zero is a valid
exponent, not a missing-value code. `target_weight` remains unchanged.

When simulation evaluation groups are written, the scalar objective in
`/evaluation/<target>/objective` receives `boxcox_lambda` and `boxcox_definition`
attributes for Box-Cox targets only. The calibration metadata does not require
full simulation time-series output.

## Apply and test

From the repository root, after the diagnostic-target patch:

```bash
git apply --check ~/Downloads/summa_configurable_boxcox_lambda_err_last.patch &&
git apply ~/Downloads/summa_configurable_boxcox_lambda_err_last.patch

git diff --check
bash utils/test/test_diagnostic_targets/run_tests.sh
bash utils/test/test_boxcox_lambda/run_tests.sh
```

This full replacement patch applies to the reviewed zero-weight diagnostic
baseline, before the original Box-Cox patch. Do not apply both full versions.

If the original Box-Cox patch is already applied, use the separately supplied
`summa_boxcox_lambda_err_last_update.patch` instead. That small migration patch
produces the same final files without reversing the old patch. Always use
`git apply --check` and inspect conflicts; do not discard unrelated edits.

Then rebuild SUMMA using the existing build procedure. Use a new calibration
output location for the integration run.

The patch changes seven production files:

* calibration/calibration_types.f90
* calibration/calibration_config.f90
* calibration/calibration_output.f90
* driver/summa_config.f90
* driver/summa_simulation.f90
* objfunc/metrics.f90
* objfunc/write_evaluation.f90

It also updates the existing diagnostic tests and writer substitute for the
new interface and adds this focused test directory. No build-system changes
or new production dependencies are needed.

## Test scope -- important

The new suite requires Python 3 and GNU Fortran, without NetCDF or TOML libraries.
It compiles the production metric module and extracts production routines
unchanged for target-scoring, parser and calibration-output tests.

It checks 182 assertions: 130 numerical/default/domain assertions, 15 scorer
assertions, 21 parser/legacy-adapter assertions and 16 metadata assertions.
The pre-existing diagnostic suite supplies another 64 assertions. The runner
also invokes `check_interfaces.py`: it checks all five affected production
interfaces, the matching test substitute, and all direct calls in the available
source/test tree for error-last ordering and keyword error arguments. It is a
focused source check, not a replacement for compiling callers or resolving
renamed/indirect procedure calls.

**These are not full SUMMA/MPI integration tests.** Observation I/O and alignment
use substitutes. Parser tests use an in-memory TOML accessor substitute (they do
not invoke the real TOML text parser). Metadata tests use a NetCDF-interface
substitute that captures calls (they do not write an actual NetCDF file). The
production evaluation writer itself is compiled and exercised in those tests.
The two example TOML files were also parsed and checked with Python `tomllib`.

An optional checker is provided for an actual calibration result:

```bash
python3 utils/test/test_boxcox_lambda/check_output.py /path/to/calibration.nc \
  --expect streamflow_nse_boxcox=0.2
```

That optional checker requires NumPy and netCDF4. It checks target metadata;
it does not recompute scores from the underlying simulation and observations.

## Optional exp50 example

The bundle contains `summa_config_template_exp50_boxcox.toml`. Relative to the
previous two-target exp50 template, its **only** configuration change is a third
zero-weight NSE/Box-Cox target with illustrative lambda 0.2. It retains ordinary
KGE at weight 1, NSE-log at weight 0, DDS, 44 evaluations, and the canopy-bottom
ratio lower bound 0.02. The original two-target template is not overwritten.

## Definition reference

The unshifted Box-Cox definition and the nonpositive-lambda domain restriction
are consistent with the SciPy special.boxcox reference:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.boxcox.html
The SUMMA implementation has no dependency on SciPy.
