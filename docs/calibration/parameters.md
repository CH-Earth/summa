# Parameter adjustments

A calibration vector contains one scalar for each sampled parameter name.
The parameterization determines how that scalar constructs a physical field.
Search bounds, search transformations, and physical-field construction are
separate operations.

Let `h` identify a selected horizontal model location and `l` a local component,
such as a soil layer. Local lengths can vary between locations. Let `p0(h,l)`
be the initialized field retained before calibration and `p(h,l)` the effective
field for a proposed trial.

| Adjustment | TOML section | Effective field | Meaning of the searched scalar |
| --- | --- | --- | --- |
| Physical value | `[calibration.parameters]` | `p(h,l) = x` | Parameter value in its physical units |
| Spatial multiplier | `[calibration.spatial_multipliers]` | `p(h,l) = m * p0(h,l)` | Coefficient scaling the initialized field |
| Parameter ratio | `[calibration.parameter_ratios]` | `p(h,l) = r * pref(h,l)` | Coefficient scaling the effective reference field |

A name may appear in only one of these three tables. Declaring a multiplier
replaces physical-value calibration of that name; it does not add a second
adjustment applied after a separately sampled physical value.

## Physical-value calibration

```toml
[calibration.parameters]
k_soil = { lower = 1.0e-7, upper = 1.0e-4, transformation = "log" }
theta_sat = { lower = 0.2, upper = 0.7 }
aquiferBaseflowExp = { lower = 1.0, upper = 4.0 }
```

If DDS proposes `k_soil = 2.0e-6`, the adapter assigns that value to every
selected element of the conductivity field. Existing spatial or vertical
variation in that field is replaced by a uniform value. Other parameters are
not made uniform merely because this one is calibrated.

This is useful for a spatially lumped basin or a controlled uniform-parameter
benchmark on a distributed domain. A uniform calibrated field does not make
the model itself spatially lumped: forcing, vegetation, topography, and other
parameters may still vary spatially.

Physical-value bounds can be omitted individually. Missing limits are obtained
from SUMMA's parameter metadata/default initialization. Explicit limits outside
the native range can generate a warning; supplying them is not a guarantee that
SUMMA's subsequent physical checks will accept every proposed value.

When both bounds are given, `lower` must not exceed `upper`. Equal bounds pin
the parameter at that value: it stays in the calibration vector and outputs,
but the search cannot move it.

A layered parameter may have no elements at some locations, such as soil
parameters in a glacier-ice domain with no soil layers. Those locations are
skipped; the parameter is rejected only when it is empty everywhere.

## Spatial multipliers

```toml
[calibration.spatial_multipliers]
k_soil = { lower = 0.1, upper = 10.0, transformation = "log" }
theta_sat = { lower = 0.8, upper = 1.2 }
```

These bounds apply to `m`, not to conductivity or porosity. Both bounds must be
specified for a multiplier. For example:

```text
initialized k_soil:  [1e-6, 2e-6, 5e-6]
proposed multiplier: 2
trial k_soil:        [2e-6, 4e-6, 1e-5]
```

A positive common multiplier preserves relative spatial and vertical variation.
A multiplier of one reproduces the initialized map, not a separately calibrated
uniform field. Multiplication cannot create variation where the initialized
field is uniform, and a zero element remains zero.

**Every candidate scales the same initialized field.** A trial with `m = 0.8`
following a trial with `m = 1.2` produces `0.8 * p0`, not `0.96 * p0`.

The initialized map is the result of SUMMA's actual initialization sequence,
including parameter tables, trial-parameter files, and applicable overrides.
Do not assume that the soil table defines the effective map if a trial-parameter
file overwrites it. Inspect the sample-zero parameter export and initialization
settings before interpreting a multiplier experiment.

### Choosing coefficient ranges

A coefficient interval is not an independent elementwise physical-bound setting.
For a positive initial field and desired physical limits `Li` and `Ui`, the
following algebra provides a useful design check:

```text
Li <= m * p0i <= Ui  for every field element i

max_i(Li / p0i) <= m <= min_i(Ui / p0i)
```

If this interval is empty, a single multiplier cannot put every element inside
that physical envelope. Reconsider the map, the desired limits, or the
parameterization; do not silently clip the map. This is a way to select TOML
bounds, not an additional range-propagation operation performed by the optimizer.
Zeros and signed fields require separate treatment.

The generic checker validates scalar bounds, finite effective fields, and the
configured ordered relationships. It does not infer every physical relationship
or silently enforce all native SUMMA limits on a multiplied field. Model-level
parameter checks still apply during initialization.

## Parameter ratios: tying one parameter to another

A ratio uses the reference parameter's **effective values for this candidate**,
not necessarily its original initialized map.

```toml
[calibration.parameters]
heightCanopyTop = { lower = 2.0, upper = 30.0 }

[calibration.parameter_ratios]
heightCanopyBottom = { lower = 0.02, upper = 0.1 }

[[parameter_dependencies.parameter_ratio]]
reference = "heightCanopyTop"
dependent = "heightCanopyBottom"
```

A proposal with canopy top `20 m` and ratio `0.05` gives canopy bottom `1 m`.
The calibration output stores `0.05` for this variable; the physical parameter
export stores the resulting height. Treating the ratio as a height would be an
incorrect replay of the trial.

Both coefficient bounds must be explicit. Ratios are not intrinsically limited
to zero through one. A relationship such as macropore conductivity being at
least soil conductivity can instead use positive ratios greater than or equal
to one, with an appropriate physical constraint if needed.

The reference can be a sampled physical value, a sampled spatial multiplier,
another ratio, or a fixed initialized parameter. It need not add a search
variable. Reference and dependent should have compatible physical meanings and
units; choosing such a relationship is a modeling decision, not something
established by the existence of the ratio interface.

References are declared in `[[parameter_dependencies.parameter_ratio]]`.
Putting `reference` inside the inline bounds object is not the documented
schema. Each dependent needs one reference; duplicate references for the same
dependent, self-references, and cycles are invalid.

### Elementwise compatibility

The default `reference_reduction = "none"` constructs elementwise ratios.
Reference and dependent must select corresponding horizontal locations and
have compatible local shapes. A layered reference does not automatically
become a scalar, nor does a scalar reference automatically broadcast into a
layered dependent through this ratio operation.

This is distinct from the local singleton broadcasting used when **checking
ordered constraints**, described in [Constraints](constraints.md).

## Chained ratios and local reference reduction

Dependencies are resolved in topological order, so a downstream ratio sees
all upstream changes. For soil-water stress, a useful explicit relationship is:

```text
effective saturation -> transpiration threshold -> wilting threshold
```

A layered saturation parameter and a scalar stress threshold need an explicit
local reduction:

```toml
[calibration.parameters]
theta_res = { lower = 0.001, upper = 0.1 }

[calibration.spatial_multipliers]
theta_sat = { lower = 0.8, upper = 1.2 }

[calibration.parameter_ratios]
critSoilTranspire = { lower = 0.5, upper = 0.95, reference_reduction = "min" }
critSoilWilting = { lower = 0.2, upper = 0.9 }

[[parameter_dependencies.parameter_ratio]]
reference = "theta_sat"
dependent = "critSoilTranspire"

[[parameter_dependencies.parameter_ratio]]
reference = "critSoilTranspire"
dependent = "critSoilWilting"

[[parameter_dependencies.ordered]]
parameters = ["theta_res", "critSoilWilting", "critSoilTranspire", "theta_sat"]
gap_fraction = 0.01
```

The effective fields are:

```text
S(h,l) = mS * S0(h,l)
T(h)   = rT * min_l S(h,l)
W(h)   = rW * T(h)
```

The minimum is taken **within each horizontal location**, before multiplying by
the ratio. It is not a minimum across the basin. The dependent must be a native
local scalar; a soil-vector parameter with one layer is not treated as a native
scalar merely because its current array length happens to be one.

For example, choose `mS = 1.1`, `rT = 0.75`, and `rW = 0.4`:

| Location | Initialized saturation | Effective saturation | Transpiration threshold | Wilting threshold |
| --- | --- | --- | --- | --- |
| A | `[0.40, 0.50]` | `[0.44, 0.55]` | `0.3300` | `0.1320` |
| B | `[0.30, 0.45]` | `[0.33, 0.495]` | `0.2475` | `0.0990` |

Three scalar coefficients produce different local threshold values. This does
not add one independent threshold coefficient per HRU.

For positive saturation and ratios strictly between zero and one, these
relationships put wilting below the transpiration threshold and that threshold
below each local saturation value. They do not establish that residual moisture
is below wilting or that every configured minimum gap is satisfied. Keep those
ordered constraints explicit.

Only `none` and `min` are implemented reference reductions in this interface.
There is no implicit root-zone average, selected-layer rule, or spatial mean.
The local-minimum construction is a deliberate parameterization choice, not a
claim that it is the best plant-water-stress model for every application.

## Remainder dependencies: closing a sum

Some parameters must add up to a fixed total, such as soil texture fractions.
A remainder dependency sets one parameter to the total minus its sources:

```toml
[calibration.parameters]
frac_sand = { lower = 0.05, upper = 0.8 }

[[parameter_dependencies.remainder]]
parameter = "frac_silt"
of = ["frac_sand", "frac_clay"]
total = 1.0
```

```text
p_silt(h,l) = total - p_sand(h,l) - p_clay(h,l)
```

`total` is optional and defaults to `1`. The remainder is evaluated at every
location and layer, after all sampled, multiplier, and ratio fields are built.
Here `frac_sand` is searched and `frac_clay` stays at its initialized values.
A source can be a sampled physical value, a multiplier, a ratio, or a fixed
parameter.

The remainder parameter adds no search variable. It must not also appear in
`[calibration.parameters]`, `[calibration.spatial_multipliers]`, or
`[calibration.parameter_ratios]`; it cannot be its own source, a source of another
remainder, or the reference of a ratio. Each parameter takes at most one remainder.
Sources must have the same horizontal locations and local lengths as the
remainder parameter.

A candidate whose sources sum to more than the total gives a negative remainder
and is rejected as infeasible (see [Constraints](constraints.md)). Remainder
values are written into SUMMA and into the trial parameter export like sampled
values.

## Search transformations

`transformation = "log"` applies to the scalar being searched. Specify bounds in
its original units:

```toml
k_soil = { lower = 1.0e-7, upper = 1.0e-4, transformation = "log" }
```

Do not supply `log(1e-7)` and `log(1e-4)` as the TOML limits. The search constructs
transformed coordinates internally and returns untransformed scalar values to
the field adapter and output writer. Log bounds must be positive. The default
search transformation is `none`.

A log-transformed multiplier is still a multiplier of the original physical
field. Parameter search transformations are unrelated to objective-series
transformations such as `obs_transform = "log"` or `"boxcox"`.

## Domain-qualified names

A quoted name such as `"k_soil@upland"` scopes an adjustment to a domain type
recognized by SUMMA. It is not an HRU-ID selector and does not mean a separate
coefficient for each selected domain.

Ratio and ordered-constraint participants must select compatible model storage
and the same domain qualifier. Do not mix GRU basin parameters with HRU/domain
parameters in one chain. Use the same qualified naming convention for the
reference and dependent. Avoid overlapping qualified and unqualified
adjustments to the same native field; this is not a supported way to define a
priority of overrides.

The current physical exporter handles domain dimensions and maps qualified names
back to their base physical parameter names. See [Outputs](outputs.md) for the
meaning and scope of these files.

## Initialization and application order

For each candidate, the field machinery restores the retained initialized
values, applies physical-value and multiplier adjustments, evaluates ratios
in dependency order, and then sets remainder parameters. Unchanged parameters keep their initialized fields. SUMMA
then updates parameter-dependent quantities and performs model-level checks.

Applying a ratio to yesterday's trial or applying calibrated parameters after
computing dependent quantities would give a different model. The adapter exists
to keep these operations explicit and in the correct order.

**Implementation entry points:** `parameter_search.f90`, `parameter_fields.f90`,
`summa_parameter_spec.f90`, `summa_parameter_fields.f90` in
`build/source/calibration/`, and `summa_setup.f90` in `build/source/driver/`.
