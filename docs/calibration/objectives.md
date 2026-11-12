# Objectives and observations

A target defines one comparison between a simulated series and an observed
series. Several targets can use the same model evaluation: there is no need to
run SUMMA again merely to calculate NSE-log alongside ordinary streamflow KGE.

## Explicit targets

```toml
[calibration]
algorithm = "dds"
n_samples = 5000
metric = "kge"
obs_transform = "none"
start_date = "2016-01-01"
end_date = "2017-01-01"
write_aligned = false

[[calibration.target]]
name = "outlet_kge"
variable = "streamflow"
spatial_unit = "reach"
spatial_id = 222
obs_path = "{home}/data/distributed/obs/"
obs_file = "CAN_07DD001_daily_flow_observations.nc"
vname_obs = "q_obs"
metric = "kge"
obs_transform = "none"
weight = 1.0
```

The reach identifier is illustrative and must match the actual routing topology
and gauging location. `spatial_id` is an identifier from the model inputs, not
an array position.

Each explicit target requires its own `obs_file`. `obs_path` and the observation
variable name can inherit defaults from `[observations]`; metric,
`obs_transform`, and `boxcox_lambda` can inherit calibration-level defaults.
A target-local exponent overrides the calibration-level exponent. An omitted
weight defaults to one.

When explicit targets exist, the legacy implicit streamflow target is **not
appended**. Include the primary KGE target explicitly when adding diagnostics.
Without an explicit target list, the observation configuration can still provide
an implicit single streamflow target; this is distinct from the removed legacy
parameter-list syntax.

Use clear, unique target names suitable for NetCDF group names. The target name
is a label; `variable`, `metric`, and spatial selection determine what is scored.

## Spatial selections and simulated variables

| Selection | Meaning |
| --- | --- |
| Omitted or `spatial_unit = "domain"` | The configured model domain; no `spatial_id` |
| `spatial_unit = "gru"` | One identified GRU |
| `spatial_unit = "hru"` | One identified HRU |
| `spatial_unit = "reach"` | One identified routed reach; streamflow/discharge targets |

`streamflow` and `discharge` are special names. A reach target uses coupled
mizuRoute discharge. A domain streamflow target uses the coupled outlet series
when routing is enabled; otherwise it sums SUMMA routed runoff times GRU area.
A GRU streamflow target uses that GRU's routed runoff and area. For a network
with several outlets or an internal gauge, use an explicit reach target rather
than assuming the default outlet represents the observation.

Other variable names are resolved against SUMMA's basin, diagnostic, prognostic,
flux, and parameter metadata. Ordinary spatial reductions use GRU or HRU area,
and HRU/domain quantities also use domain-area fractions. Zero-area domains are
excluded when a positive domain-area total is available. The code has fallback
weights for nonpositive areas; valid area metadata are therefore essential to
interpreting an aggregate physically.

**Naming a layered variable does not request a vertical profile or a root-zone
average.** The current ordinary HRU/domain value accessor reads its first local
component, `dat(1)`, then performs the spatial aggregation. Use an appropriate
scalar diagnostic for integrated or depth-specific observations, or extend the
collector explicitly. Do not assume arbitrary vertical observation operators.

The same variable at the same spatial unit is collected once and reused by
multiple targets. The same variable at different locations requires separate
simulated series.

## Observation files and time conventions

The observation reader expects a NetCDF `time` dimension and coordinate,
a one-dimensional observed series named by `vname_obs`, and units for time and
the observed variable. Missing-value metadata are handled by the reader.
An optional `time_bnds` variable supplies the observation intervals.

When bounds are absent, observations are interpreted as covering the interval
since the preceding timestamp. The first interval uses the spacing to the next
observation. Thus, timestamps are effectively period-ending in this fallback
convention. A single timestamp cannot establish a nonzero interval by itself.
Provide explicit bounds when periods are irregular or the observation product
uses a different reporting convention.

For ordinary targets, simulated values are averaged over each observation's
interval using model timestamps in `(interval_start, interval_end]`. The aligned
record retains the number of contributing simulation steps. This is not a
universal integration or interpolation scheme for arbitrary observations.
Check coverage at the first and last intervals and across missing periods;
a contributing-step count does not, by itself, prove full interval coverage.

The calibration date filter uses observation timestamps in the half-open window:

```text
start_date <= observation_time < end_date
```

The example above selects observation timestamps from the beginning of 2016 up
to, but not including, the beginning of 2017. For daily period-ending data,
decide explicitly whether a timestamp on January 1 represents the preceding
day and choose the window accordingly. Do not assume that an `end_date` names
an inclusive final day.

Time coordinates are converted using their supported units/reference time.
Observed and simulated value units must be compatible. The optional
`obs_units` setting overrides the interpreted units label; it is not a general
numerical unit-conversion facility. For example, relabeling cubic feet per
second as cubic metres per second would not convert the values. The current
reader still expects the observation variable's units attribute to be present.

## Available metrics

| Metric | Interpretation | Preferred direction |
| --- | --- | --- |
| `kge` | Kling-Gupta efficiency using correlation, standard-deviation ratio, and mean ratio | Maximize |
| `kgep` | Modified KGE using coefficient-of-variation ratio | Maximize |
| `nse` | Nash-Sutcliffe efficiency | Maximize |
| `mae` | Mean absolute error | Minimize |
| `rmse` | Root mean square error | Minimize |

For paired transformed observations `o` and simulations `s`, the efficiency
components are:

```text
beta  = mean(s) / mean(o)
alpha = sd(s) / sd(o)
gamma = [sd(s) / mean(s)] / [sd(o) / mean(o)]

KGE  = 1 - sqrt((cor(o,s)-1)^2 + (alpha-1)^2 + (beta-1)^2)
KGE' = 1 - sqrt((cor(o,s)-1)^2 + (gamma-1)^2 + (beta-1)^2)
NSE  = 1 - sum((s-o)^2) / sum((o-mean(o))^2)
```

Metrics operate on usable paired values. Insufficient pairs, required zero
variances or means, invalid transform domains, and non-finite results can make a
metric undefined. Very poor but finite scores are not clipped to an artificial
minimum. Undefined scores are handled differently for active targets and
zero-weight diagnostics, as described below.

## Transforming observed and simulated values

Despite its name, `obs_transform` transforms **both** observations and simulations
before computing the metric. It is separate from parameter search transformations.

| Setting | Transformation |
| --- | --- |
| `"none"` | No value transformation |
| `"log"` | `log(Q + epsilon)`, with epsilon derived from valid observations |
| `"boxcox"` | Unshifted Box-Cox transformation with configured lambda |
| Numeric string, for example `"0.5"` | Power transformation with that exponent |

The logarithmic offset is `mean(valid observations) / 100`, applied identically
to observation and simulation values. It is data-dependent, not a global
constant. Negative power exponents also use the observation-derived offset in
the implemented power path. The transformed arguments must satisfy the relevant
domain checks. Literal names such as `"sqrt"` are not substitutes for the numeric
power string in this interface.

### Box-Cox exponent

```text
B(Q; lambda) = (Q^lambda - 1) / lambda    when lambda != 0
B(Q; 0)      = log(Q)
```

The default is `0.25`; lambda is fixed for evaluation, not sampled by DDS. The
implementation preserves the historical arithmetic at the default and handles
the near-zero limit without unstable subtraction of nearly equal quantities.

For positive lambda, zero flow is allowed. Nonpositive lambda requires strictly
positive flow. Negative flow is invalid for this Box-Cox implementation. No
epsilon is added: `boxcox_lambda = 0.0` is **not** equivalent to the shifted
`obs_transform = "log"` option.

## DDS diagnostics without changing optimization

Append explicit diagnostic targets while retaining the explicit primary target:

```toml
[[calibration.target]]
name = "outlet_nse_log"
variable = "streamflow"
spatial_unit = "reach"
spatial_id = 222
obs_path = "{home}/data/distributed/obs/"
obs_file = "CAN_07DD001_daily_flow_observations.nc"
vname_obs = "q_obs"
metric = "nse"
obs_transform = "log"
weight = 0.0

[[calibration.target]]
name = "outlet_nse_boxcox_025"
variable = "streamflow"
spatial_unit = "reach"
spatial_id = 222
obs_path = "{home}/data/distributed/obs/"
obs_file = "CAN_07DD001_daily_flow_observations.nc"
vname_obs = "q_obs"
metric = "nse"
obs_transform = "boxcox"
boxcox_lambda = 0.25
weight = 0.0
```

DDS skips zero-weight targets before accessing their metric values during
scalarization. An undefined diagnostic is therefore allowed to remain missing
without corrupting a valid optimization score. An undefined active metric marks
the trial as failed. Observation, configuration, alignment, and output errors
are not made harmless simply by assigning zero weight to the affected target.

Zero-weight diagnostic targets are currently rejected for NSGA-II. Positive
weight magnitudes do not control NSGA-II Pareto selection.

All target scores are written to the calibration file even with
`write_timeseries = false`. `write_aligned` controls optional detailed evaluation
series, not whether the calibration objectives are calculated.

## Accumulation and storage anomalies

Use `accumulate = true` when a simulated rate must first be integrated into a
cumulative series. This occurs before observation alignment. The implemented
path includes the storage-change use case, converting the supported water-mass
rate units into cumulative water-equivalent depth; it is not a declaration that
every arbitrary flux unit can be integrated and converted automatically.

For ordinary anomaly targets, supply both `baseline_start` and `baseline_end`.
The driver aligns over the union of evaluation and baseline periods, identifies
matching finite observation/simulation periods, subtracts a separate mean from
each series using that common support, then retains only the evaluation window.
The two means need not be numerically equal; their temporal support must match.

For example, a storage-change target can use `variable = "basin__StorageChange"`,
`accumulate = true`, and `metric = "rmse"` against a water-storage anomaly
observation series. Confirm the variable meaning, units, spatial footprint, and
baseline period for the actual product.

`balance = "seasonal"` or `"annual"` activates the specialized balance-between-
extremes calculation. `balance_window` controls the search window around its
expected extrema. This is a specialized observation operator, not general
calendar-period differencing; inspect `balance_between_extremes` before applying
it to a new balance dataset. Its baseline processing follows the specialized
balance path rather than the ordinary matched-anomaly path.

## Low-flow subsets and flow-duration comparisons

The checked merged implementation additionally accepts `low_flow` and
`flow_duration` in a target. These are series-selection operations, not metric
names or optimizer settings:

```toml
# Additional fields within an existing [[calibration.target]] block:
low_flow = 0.2
flow_duration = false
```

`low_flow` defaults to `1.0` and must be greater than zero and no greater than
one. With `flow_duration = false`, the code uses the lower fraction of the
observed finite aligned values and retains their corresponding simulated
values. Ties at the observation cutoff can retain more than the nominal
fraction. Temporal correspondence is preserved.

With `flow_duration = true`, the valid aligned observation and simulation values
are sorted independently and the requested lower fraction of each sorted series
is compared. At `low_flow = 1.0` this compares their full sorted distributions.
Temporal correspondence is intentionally lost. This is different from choosing
low-flow dates using observations and comparing simulated flow on those dates.

Selection occurs after alignment/baseline processing and before the metric's
value transformation. In particular, a shifted logarithm then derives its offset
from the retained observed values. These settings are not fully represented in
the calibration target metadata written by the current output module: retain
the TOML when using them.

## Scientific interpretation

A diagnostic target does not independently calibrate another parameter set.
For a KGE-focused DDS experiment, select the sample using KGE and report all
other metrics at that same sample. Optimizing each displayed diagnostic
separately would answer a different question.

Similarly, improving an outlet metric is not sufficient evidence that every
internal flux or storage improved. Use independent periods, additional gauges,
and process-relevant observations to evaluate the selected simulations where
those data are available.

**Implementation entry points:** `calibration_config.f90` and
`calibration_types.f90` in `build/source/calibration/`; `summa_simulation.f90`
in `build/source/driver/`; and `read_observations.f90`, `time_align.f90`,
`metrics.f90`, `series_transform.f90`, and `simulated_series.f90` in
`build/source/objfunc/`.
