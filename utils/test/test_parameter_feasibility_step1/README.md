# Step 1: nested physical fields and generic feasibility

## Status and scope

Step 1 is implemented and tested independently of the production optimizers.
The candidate vector remains generic; physical fields now use:

```fortran
field(iParam)%initial%hspatial(iSpatial)%dat(iLocal)
field(iParam)%effective%hspatial(iSpatial)%dat(iLocal)
```

There are no flattened field arrays, block IDs, offsets, expanded broadcasting
arrays, or feasible-interval solver in the new `src/` modules. `dat(:)` is a
local component vector, not necessarily soil depth.

**Do not copy these four files over the production modules yet.** Their field
representation and construction interfaces are deliberately different. Step 2
will connect the existing sampler rejection loops to the checker and update
production callers. This isolated package establishes the new core first,
without breaking the currently working calibration branch.

The optional additive patch installs this experiment under
`utils/test/test_parameter_feasibility_step1/`; it does not modify production
sources, CMake, DDS, NSGA-II or MPI dispatch.

## Sources and provenance

The existing parameter specification/search types and mathematical semantics
come from the uploaded `summa_calibration_source.tar(1).gz`, the later source
files/diffs, and the earlier `summa_feasibility_design` prototype. The code here
implements the subsequent agreed nested representation. This is not a full
current SUMMA checkout.

`support/` contains UNMODIFIED uploaded dependencies for standalone tests:
`nr_type.f90`, `parameter_search.f90`, `var_lookup.f90`, `calibration_types.f90`,
and `data_types.f90`. Their original headers are retained. Their hashes and the
new files' hashes appear in `manifest.json`.

The production source change is the staged code in `src/`, not the old algorithms
in `support/`. The test-only helper `build_search` calls the existing search
initializer on a COPY of the specification without its ordered chains. The
original specification retains every chain. This bypasses the old scalar bound
propagation while testing the new pure checker, without editing the optimizer.
This bridge is for Step 1 only; it is not the proposed production integration.

## Run

After extracting the bundle:

```sh
cd summa_feasibility_step1
./run_tests.sh
```

Requires GNU Fortran and a POSIX shell. There are no NetCDF or MPI dependencies
for these tests. `FC` and `FFLAGS` can be overridden:

```sh
FFLAGS='-O2 -g -fcheck=all -fbacktrace -Wall -Wextra' ./run_tests.sh
```

Do not enable `-ffast-math` or `-ffinite-math-only`: the checker deliberately
handles NaN/Inf candidates. Overflow tests use IEEE nontrapping behavior; do
not add floating-point overflow traps to this particular test suite.

The runner builds in a temporary directory and leaves no `.o` or `.mod` files
in the source tree. In the standalone bundle it uses `support/`; when installed
by the additive patch (without `support/`), it uses the checkout's shared source
files. An actual full SUMMA build is not part of this runner.

## Four staged modules

| Module | Responsibility |
| --- | --- |
| `parameter_fields` | Nested types, setup validation, resolved ratio references/evaluation order, construction/reset of native-sized fields |
| `parameter_constraints` | Configured scalar bounds, finiteness, physical gaps, all ordered comparisons and local singleton broadcasting |
| `model_calibration` | Generic context, setup entry point and opaque `is_feasible` callback |
| `summa_parameter_fields` | SUMMA extraction/correspondence checks and accepted-field write-back |

The only extra construction metadata is `reference(:)` and `order(:)`, resolved
once in a `parameter_field_plan`. These arrays index PARAMETER DEPENDENCIES;
they do not represent spatial data or replace the discarded block machinery.

All routines with error outputs retain `err,message` last.

## Public generic interface

Populate the context's specification, raw configured search mapping/bounds and
initialized fields, then initialize once:

```fortran
call initialize_model_constraints(context,err,message)
```

This prepares the field construction order, allocates effective arrays, validates
local shapes and all constraint definitions, and checks entirely fixed chains.
It does not evaluate a lower-bound candidate or refine parameter bounds.

Every complete untransformed candidate uses the same call:

```fortran
call is_feasible(context,x,feasible,err,message)
```

- `err=0`, `feasible=.true.`: accept.
- `err=0`, `feasible=.false.`: ordinary candidate rejection.
- `err/=0`: configuration/API/checking failure; do not resample this error.

The optimizer need not interpret the message or know which constraint failed.
A non-finite constructed candidate field is a rejection. A non-finite INITIAL
field is a setup error. No retries or candidate generation occur inside the
checker.

## Representation and correspondence

```fortran
type :: parameter_values
  real(rkind), allocatable :: dat(:)
end type

type :: spatial_parameter_field
  type(parameter_values), allocatable :: hspatial(:)
end type

type :: parameter_field
  type(spatial_parameter_field) :: initial
  type(spatial_parameter_field) :: effective
end type
```

Example:

```text
theta_sat:          hspatial(1)%dat = [0.40, 0.42, 0.44]
                    hspatial(2)%dat = [0.38, 0.41]
critSoilTranspire:  hspatial(1)%dat = [0.30]
                    hspatial(2)%dat = [0.35]
```

For an adjacent pair, at a corresponding location with local lengths `na` and
`nb`, compare `a%dat(min(k,na))` and `b%dat(min(k,nb))`, for
`k=1:max(na,nb)`. Setup permits only equal local lengths or local 1-to-N.
There is NO horizontal singleton broadcasting. Empty arrays are setup errors.

The generic structure does not encode physical identities or axis names.
An adapter must establish that corresponding `hspatial(i)` entries represent
the same location and corresponding nonscalar `dat(k)` entries the same local
support. Equal array sizes cannot establish that.

The supplied SUMMA adapter enforces consistent GRU-versus-HRU storage and domain
qualifiers for linked parameters. It uses existing `mpar_meta%varType` to reject
incompatible nonscalar local axes (such as soil versus spectral vectors), even
when lengths coincide. A true `scalarv` parameter may be compared with a vector.
It traverses the selected GRU/HRU/domain entries in one deterministic order.
No such model-specific interpretation enters the generic constraints module.

Ratio CONSTRUCTION retains the existing same-native-local-length requirement;
this step adds no new scalar-to-layer ratio construction rule. Local scalar
broadcasting is for ORDERED COMPARISONS only.

## Preserved physical semantics

For every location/local element:

```text
physical:    P = x
multiplier:  P = x * P_initial
ratio:       P = x * P_reference_effective
fixed:       P = P_initial
```

References may be fixed or sampled and may themselves be ratios. The resolved
evaluation order handles arbitrary registry order; cycles are setup errors.
Every candidate resets effective values from initial fields. No accumulation
of multipliers or ratios occurs.

For each adjacent ordered pair:

```text
P_b - P_a >= gap - tolerance
```

The whole-field gap definition is unchanged:

```text
gap = gap_fraction * (maximum_attainable_last - minimum_attainable_first)
```

Extrema include ALL horizontal locations and ALL local components. They are not
computed separately for each HRU. Fixed endpoints use their actual initialized
extrema; sampled physical endpoints use configured scalar bounds; multiplier
endpoints use bounds times initial-field extrema; ratio endpoints use bounds
times CURRENT effective reference-field extrema. Constant gaps are cached;
ratio-endpoint gaps are recomputed for the candidate.

A zero gap permits equality within tolerance. The magnitude-scaled tolerance
is retained from the earlier physical-field validator. Fixed-parameter nominal
bounds do not replace their actual initialized physical fields.

Search bounds are NOT tightened. Setup rejects a mismatch between configured
bounds and the supplied search intervals to avoid silently retaining the old
scalar-bound propagation in the new pure-rejection design.

## Ownership and application

The specification, raw search bounds/mapping, initialized fields, shapes and
spatial/local ordering are immutable after setup. Reinitialize after editing
any of them. The checker mutates ONLY its context-owned effective scratch.
Use one context per concurrent checker/case. Do not share it concurrently
between threads. Inspect effective fields for application only after acceptance.

The adapter's `apply_summa_calibration_adapter` applies an ALREADY ACCEPTED
vector without repeating calibration inequalities. Its caller must first check
that vector with the same initialized context. It still checks finite effective
values and every destination shape BEFORE performing any model write. Scalar
native arrays remain scalar. It contains no validation-only/prevalidated flags.

A model's independent physical/numerical checks remain necessary. Calibration
feasibility is not a guarantee of a successful hydrological simulation.

## Tests actually executed

With GNU Fortran 14.2.0, both `-O0` and `-O2`, plus runtime bounds/allocation checks:

- Generic nested core: 4,190 assertions passed.
- SUMMA nested adapter: 61 assertions passed.
- Separate comparison with the earlier block prototype: 5,683 candidate
  decisions matched exactly; constructed-field sums matched within 1e-14
  relative/absolute tolerance (nested reductions change addition order).

These tests include synthetic soil and spectral examples, unequal local lengths,
all-fixed configurations, scalar/vector comparisons in both directions,
whole-field and dynamic ratio gaps, reference chains/cycles, signed/zero scales,
NaN/Inf/overflow, malformed shapes/mappings, and preservation of native arrays.
The adapter uses REAL SUMMA data/mapping/metadata types; only globalData is
replaced with a small test holder. No model equations or MPI execution are stubbed
and then claimed as tested: they are simply outside these tests.

The compiler warnings are real-equality checks in unchanged dependencies and an
intentional exact configured-bound identity check in the new field initializer.
Logs are in `test_results/` in the standalone bundle.

Optional parity rerun, with the earlier prototype extracted separately:

```sh
./compare_reference.sh /path/to/summa_feasibility_design
```

This optional comparison also requires Python 3.

NOT executed: full SUMMA build, DDS/NSGA-II callback integration, MPI dispatch,
or Athabasca hydrological runs. Those are Step 2/integration work.

## Next integration boundary

The future sampler acceptance call is `is_feasible`; existing generator retry
loops will own rejection. Step 2 will replace their old scalar-only checks,
remove the newer outer dispatcher loops, and update production application
callers to the nested representation. Do not combine this isolated core with
the current production samplers by overwriting source files.

The new `src/` modules are measured separately in `source_size.json`. A staged
side-by-side test package adds files; it is not yet a net production-repository
LOC reduction. That reduction is realized when superseded production paths are
removed during integration.
