# Outputs and diagnostics

Calibration output serves three different purposes: recording the search,
exporting selected physical parameter fields, and retaining enough provenance
to identify a run. These are related but are not interchangeable restart formats.

## Filenames and isolation

The normal calibration file prefix contains the case name and independent run
index, with any configured output suffix included:

```text
<output_path>/<case_name>_run0001[_suffix]_calibration.nc
<output_path>/<case_name>_run0001[_suffix]_sample0_trialParams.nc
<output_path>/<case_name>_run0001[_suffix]_sample<ID>_trialParams.nc
```

The experiment tag is not automatically part of every filename. Separate
experiment/job output directories are therefore necessary when two jobs may
produce the same case/run names. Files are created for a new search, not opened
as an optimizer continuation.

## Calibration NetCDF schema

| Variable or metadata | Interpretation |
| --- | --- |
| `sample` | Proposed trial identifier |
| One variable per sampled parameter | Untransformed scalar calibration variable for each trial |
| `objective` | One metric value per target per trial |
| `target_name` | Target labels in objective order |
| `target_weight` | Configured weights; zero means DDS diagnostic only |
| `target_boxcox_lambda` | Effective exponent for Box-Cox targets; fill for other transforms |
| `trial_failed` | Failed versus successful evaluation |
| `worker_rank` | Worker rank inside the case's instance communicator |
| `dispatch_time`, `completion_time` | Eight-component Fortran date/time records |

The objective array has both sample and target dimensions. Fortran defines the
array with target as its leading index; other readers may display dimensions
in another order. Inspect dimension names rather than flattening the array or
assuming the first column is the desired metric.

Metric names, transformations, simulated-variable names, and preferred directions
are recorded in target order on the objective variable. Read these with
`target_name` and weights to distinguish, for example, ordinary KGE from a
zero-weight NSE-log diagnostic.

The `sample` dimension is allocated to the requested budget before completion.
Its length is not evidence that every trial ran. Use valid completion records
and failure status to distinguish completed, failed, and unwritten samples.
The eight timestamp components are:

```text
year, month, day, UTC_offset_minutes, hour, minute, second, millisecond
```

Trial IDs describe proposal order, not asynchronous completion order. Use the
recorded times when constructing wall-clock convergence curves. Worker rank one
in two different case groups does not identify the same global MPI process.

### Parameter values are not always physical values

Each sampled parameter variable includes `param_type`, bounds, transformation,
and where applicable `reference_name` and `reference_reduction`.
The global `parameter_space` is `untransformed_calibration_variables`.

For a physical parameter the stored scalar is a physical value. For a spatial
multiplier it is a coefficient; for a ratio it is a tied-parameter coefficient.
A log search transformation is already inverted in the stored scalar. Do not
exponentiate it again, and do not treat every parameter-named variable as an
HRU-level physical parameter.

Ordered-constraint names and gap fractions are also recorded. The historical
scalar `trial_value` representation is not a substitute for the retained
spatial initialization fields.

## Selecting the result to analyze

For DDS, reconstruct the configured direction-aware weighted scalar from valid
active objectives and select its best completed, successful sample. With one
active ordinary KGE target, this reduces to maximizing that KGE column.

Report zero-weight diagnostics at **the same selected sample**. Selecting the
best KGE sample, best NSE-log sample, and best Box-Cox sample independently would
combine three different parameter sets and overstate the performance of any one
calibrated model.

An undefined zero-weight diagnostic is stored as the objective's fill value
while `trial_failed` remains zero. Do not discard a successful KGE trial merely
because an unrelated diagnostic is missing. A failed trial has its objective
vector filled and cannot be a valid winner.

For comparisons across configurations, match basin identities and use consistent
budgets and repeat-selection rules. A best-of-two result has a different
selection opportunity from an individual run. Incomplete cases should be
reported as incomplete, not silently interpreted as poor completed scores.

## NSGA-II output

NSGA-II additionally records:

| Variable | Meaning |
| --- | --- |
| `birth_generation` | Generation in which a trial was proposed |
| `population` | Sample identifiers retained after each generation |
| `population_rank` | Nondomination rank for retained population members |
| `population_crowding` | Crowding distances for retained members |
| `pareto_front` | Flags identifying nondominated trials across evaluated trials |

Orient each target using its metric direction when checking dominance. Do not
use weight magnitudes as Pareto trade-off coefficients. Select a compromise from
the Pareto alternatives according to an explicit application criterion, rather
than assuming the file contains a uniquely best parameter set.

## Physical parameter exports

`sample0_trialParams.nc` records the initialized baseline before calibration; it
is not an evaluated trial and does not consume a sample. Selected trial files
contain physical values after applying the parameterization. DDS selects its
elite; NSGA-II exports selected nondominated trials.

The exporter combines parameters present in the original trial-parameter file
with sampled parameters. It is not a dump of every internal state or every
possible model parameter. Parameters outside that export still depend on the
same model tables/defaults during replay.

The checked merged exporter supports multiple domains within HRUs and
maps domain-qualified calibration names to their base physical variable names.
It writes domain and supported soil-layer dimensions and missing entries where
required by the rectangular NetCDF layout. This replaces the older one-domain
export restriction. Always preserve the corresponding HRU/GRU/domain mapping
and parameter metadata when moving the exported file to another configuration.

A standalone replay requires the exported physical parameter file **and** the
appropriate initial state, simulation period, decisions, tables, forcing,
attributes, and routing configuration. Reusing only the scalar coefficients from
the calibration file through physical `--param` overrides is not a valid replay
for multiplier or ratio parameters.

## Optional evaluation groups

When model time-series output is enabled, objective diagnostics can also be
written to `/evaluation/<target_name>` in the model output. `write_aligned =
true` additionally requests aligned `time`, `time_bnds`, `obs`, `sim`, and
`n_sim_steps` data there. A compatible timestep output file must be available.

These groups do not replace the calibration NetCDF. With
`write_timeseries = false`, scores are still calculated and written to the
calibration file, but these detailed model-output groups are not written.
Aligned values describe the alignment stage; they are not necessarily already
value-transformed or reduced to a low-flow subset used by the final score.
Retain target settings to reproduce those later operations.

## Minimal provenance

The calibration file records the resolved decisions-file path, case-leader
hostname, compiled Git hash, independent run index, and cases per node:

```text
model_decisions_file
compute_hostname
git_hash
run_index
cases_per_node
```

Existing attributes include case name, worker count, node index, global case
group, and algorithm. The Git value comes from the executable's build stamp,
not a runtime query of an arbitrary source checkout. Reconfigure/rebuild when
that stamp needs updating.

A pathname is not a content fingerprint. These attributes do not embed all
inputs, an expanded TOML, a reliable dirty-tree record, or an experiment
identifier. Archive the manifest/template and actual input versions alongside
production results, especially for settings not fully represented in output
metadata, such as low-flow selection options.

## Failures and persistent output

Recoverable trial failures are recorded and excluded from elite selection. The
worker also writes local failed-trial information, including scalar calibration
values and the error message, with configuration and initial-state material for
investigation. These records may contain paths from node-local storage.

Each failed-trial directory also holds a `rerun.sh`, which replays the trial
with the serial driver by passing the trial's values as `--param` overrides.
Replay works when the case was launched from a single TOML file and every
calibrated parameter is a physical value; domain-qualified names are passed
through, so the physical fields are reproduced. Spatial multipliers and
parameter ratios cannot be expressed as `--param` values, and a manifest case is
built from a template, so in those cases `rerun.sh` prints the reason and exits.
Rebuild the physical fields and initialization from `trial.json` by hand.

The current evaluator also classifies errors from trial initialization broadly
as failed trials. Consequently, repeated failures with the same initialization
message should be investigated as possible configuration or infrastructure
problems, not assumed to be harmless bad parameter samples. Target-scoring
input/alignment errors are distinct from mathematically undefined metrics.

`--persistent-output` copies the calibration file and selected physical parameter
exports after a case finishes normally. It does **not** automatically copy the
sample-zero file, common spin-up restart, failed-trial directories, or all logs.
Temporary case logs are removed after normal cleanup. Preserve additional files
explicitly when required for diagnosis or exact replay.

Persistence at case completion is not a checkpoint system. If a job stops before
a case finishes, its useful partial results can remain only on ephemeral storage.
Completed cases in a manifest job can already have persistent copies while the
last active cases do not.
