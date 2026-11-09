# Constraints and feasibility

The optimizer proposes a vector of scalar calibration variables. Feasibility is
assessed on the **effective physical fields constructed from that vector**, not
just on the scalars as written in the TOML file.

This matters when one scalar controls many layers or locations, when references
are tied through ratios, or when a sampled parameter is constrained against a
fixed spatial field.

## Ordered parameter chains

Declare the parameters from lower to higher:

```toml
[[parameter_dependencies.ordered]]
parameters = ["theta_res", "critSoilWilting", "critSoilTranspire", "theta_sat"]
gap_fraction = 0.025
```

Every adjacent pair must satisfy:

```text
p[j+1](h,l) - p[j](h,l) >= gap
```

A small floating-point tolerance is used in the comparison. With
`gap_fraction = 0`, the condition is nondecreasing order, allowing equality within
tolerance; it is not a strict less-than constraint.

An ordered declaration does **not** automatically make all its members sampled.
A member absent from the three calibration parameter tables remains fixed at its
initialized physical values, while still participating in the constraint. This
is a major distinction from the removed `sample_gaps` mechanism.

## What `gap_fraction` means

The minimum separation is based on the envelope from the first member's lower
physical limit to the last member's upper physical limit:

```text
gap = gap_fraction * (U_last - L_first)
```

That same gap applies to every adjacent pair in the chain and every selected
location. For example, if the first lower limit is `0.001`, the last upper limit
is `0.7`, and `gap_fraction = 0.025`, then the minimum adjacent gap is `0.017475`.
It is not 2.5 percent of the current local saturation and not 2.5 percent of
each individual adjacent interval.

The endpoints are obtained according to the parameterization:

| Endpoint type | Physical envelope used for the gap |
| --- | --- |
| Sampled physical value | Its configured scalar physical bounds |
| Fixed parameter | Minimum and maximum of its initialized field |
| Spatial multiplier | Envelope of coefficient bounds multiplied by initialized field extrema |
| Ratio | Envelope of coefficient bounds multiplied by effective reference extrema |
| Ratio with local `min` reduction | Envelope using extrema across the local reference minima |

For products, the implementation considers all lower/upper coefficient and
minimum/maximum field combinations. A ratio endpoint can make the gap depend
on the current candidate because its reference field is effective, not fixed.
This calculation defines the constraint; it does not rewrite the search bounds.

Choose compatible units for chain members. A gap has the units of the physical
values being compared. Ordering unrelated quantities is not meaningful simply
because their numerical values can be compared.

## Spatial and local compatibility

The SUMMA adapter establishes correspondence between the selected locations.
Chain members must belong to compatible model storage and domain selection.
A chain cannot mix a GRU basin parameter with an HRU/domain parameter, or compare
upland-only values against a differently selected horizontal field.

At one horizontal location, an ordered comparison can compare a local scalar
against every element of a compatible local vector. Two nonscalar vectors must
have corresponding local axes and compatible lengths. Equal array lengths alone
do not make unrelated quantities, such as soil layers and spectral bands,
compatible.

There is no automatic broadcasting across horizontal locations. Likewise, local
singleton broadcasting for a constraint does not relax the elementwise shape
requirements for constructing a parameter ratio.

## The feasibility sequence

Before evaluating a candidate, the generic path:

1. Checks the scalar vector length, finite values, and configured search bounds.
2. Restores initialized fields and constructs physical, multiplier, and ratio
   adjustments in dependency order.
3. Rejects invalid or non-finite constructed values.
4. Checks every candidate-dependent ordered chain over its selected physical
   locations and local components.

Fixed-only ordered chains are checked at initialization. An invalid fixed-only
chain cannot be repaired by drawing another candidate, so it is reported as a
setup problem rather than allowing a futile search.

The shared callback contract distinguishes an ordinary rejection from a failure
of the configuration or checker:

```text
feasible = true,  err = 0 : candidate may be evaluated
feasible = false, err = 0 : reject and propose another candidate
err != 0                 : configuration/checker error; do not silently retry
```

DDS/random proposals and NSGA-II offspring use this same predicate. Each search
algorithm retains its own finite retry loop. A very restrictive feasible region
can therefore make proposal generation expensive or exhaust the retry limit.
Reaching that limit is not a reason to disable physical constraints without
examining their compatibility.

## What is not checked generically

The generic feasibility machinery does not encode every SUMMA physical rule.
It checks the registered bounds, finite constructed fields, spatial compatibility,
and declared ordering. SUMMA's parameter checks and model physics can still
reject a candidate during initialization or integration.

In particular, coefficient bounds do not create an automatic elementwise clamp
on physical fields. A porosity multiplier can pass its scalar bounds while
producing a physical porosity that a later model check rejects. Select sensible
coefficient intervals and examine failed-trial records.

## Rejected proposals versus failed trials

An infeasible proposal rejected before dispatch is not a model evaluation and
is not written as a completed calibration trial. It does not consume the
`n_samples` model-evaluation budget.

A dispatched candidate that fails during model initialization, integration, or
active-metric evaluation is a failed trial. It consumes a trial index, is recorded
with `trial_failed`, and cannot become the DDS elite. This distinction matters
when comparing useful evaluations, wall time, and failure rates.

## Removed behavior

There is no gap-coordinate search in this interface. `sample_gaps` is rejected,
not ignored. There is also no automatic scalar bound-propagation pass that
changes configured search limits to make every proposal ordered by construction.
The optimizer searches the configured coordinates and accepts candidates through
physical-field feasibility checks.

Switching from the previous bound/gap machinery changes the search problem seen
by the optimizer. Identical seeds and budgets should not be interpreted as a
promise of identical trajectories across that migration.

**Implementation entry points:** `parameter_constraints.f90`,
`parameter_fields.f90`, `model_calibration.f90`, and
`summa_parameter_fields.f90` in `build/source/calibration/`.
