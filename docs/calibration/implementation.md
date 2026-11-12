# Implementation and extension

The refactor separates generic search/field operations from SUMMA storage and
execution. This page describes the current division of responsibility, not a
claim that the whole calibration workflow is already a standalone library.

## Module map

Paths below are relative to `build/source/`.

| Module/file | Responsibility |
| --- | --- |
| `calibration/calibration_types.f90` | Calibration and target configuration types |
| `calibration/calibration_config.f90` | Calibration TOML parsing and basic settings |
| `calibration/parameter_search.f90` | Parameter specification, transformed coordinates, random and DDS proposals |
| `calibration/parameter_fields.f90` | Nested initialized/effective fields and ratio dependency order |
| `calibration/parameter_constraints.f90` | Ordered physical-field checks and endpoint gap envelopes |
| `calibration/model_calibration.f90` | Generic context and optimizer-facing feasibility callback |
| `calibration/nsga2.f90` | Evolutionary operators, nondomination ranking, and population selection |
| `calibration/summa_parameter_spec.f90` | Resolve SUMMA names, metadata, bounds, qualifiers, and fixed references |
| `calibration/summa_parameter_fields.f90` | Map native SUMMA storage to fields and apply accepted samples |
| `calibration/summa_parameter_sampling.f90` | Coordinator/worker search execution and MPI trial messages |
| `calibration/calibration_parallel.f90` | Node and instance communicator setup |
| `calibration/dynamic_case_dispatch.f90` | Generic global case queue and callback |
| `calibration/calibration_output.f90` | Trial/generation records and failed-trial diagnostics |
| `calibration/summa_parameter_export.f90` | Export initialized/selected native physical parameters |
| `driver/summa_driver_opt.f90` | Program entry, configuration mode, MPI lifecycle |
| `driver/summa_calibration_case.f90` | Execute one assigned basin/run case |
| `driver/summa_spinup.f90` | Common spin-up and initial parameter fields |
| `driver/summa_simulation.f90` | Simulation lifecycle, target scoring, and metric orientation |
| `objfunc/` modules | Observations, simulated-series collection, alignment, transformations, metrics, and evaluation output |

The generic files avoid knowledge of GRU/HRU storage. The scheduler avoids
SUMMA-specific case internals. Model evaluation, parameter export, and parts of
trial orchestration/output remain model-specific integration code.

## Case lifecycle

The execution sequence is:

```text
summa_driver_opt
  initialize command/configuration mode
  initialize world/node/instance communicators
  direct case, or generic dynamic case queue
    execute_summa_case(context, case_index, err, message)
      run_case
        resolve basin name and independent run index
        read template/configuration and expand paths
        check search settings
        prepare per-case filenames and logging
        spin up baseline; retain mpar_init and bpar_init
        initialize parameter specification, search coordinates, and coordinator seed
        initialize field adapter and constraint plans
        export sample-zero physical parameters
        create calibration NetCDF and write minimal provenance
        coordinator dispatches; workers evaluate
        export selected physical parameter fields
        close/synchronize output and clean up
        persist selected results
  finalize MPI
```

Some setup work is performed on every rank; coordinator-only file writing is
explicit. The instance communicator carries candidate/result messages. The
single-rank domain communicator is used for each model evaluation.

## Explicit case context

`summa_case_context` carries a pointer to the driver configuration, the instance
parallel context, node index, and global case-group index. A module-level
`execute_summa_case` validates the context type and configuration association
before calling `run_case`.

The generic dispatcher accepts this context as `class(*)` and a callback; it
does not import SUMMA's configuration or model types. This removes the previous
internal-procedure callback that captured host variables and makes the boundary
between scheduling and model execution explicit.

The context is rank-local state. Passing a Fortran object to a local callback is
not equivalent to serializing that object across MPI. The caller must maintain
the lifetime of its configuration target throughout dispatch.

## Scalar registry and nested fields

`parameter_spec` contains entries for sampled variables and any additional fixed
members needed by references/constraints. `parameter_search_info` maps registry
indices to the compact sampled vector and retains both original and transformed
bounds.

The generic physical representation uses horizontal entries with local arrays,
conceptually `hspatial(h)%dat(:)`. It does not require a uniform layer count or a
flat global scalar per parameter. `parameter_field` retains separate initialized
and effective values. The field plan records dependency order and local
reference-reduction choices.

`summa_parameter_map` identifies native model or basin parameter storage and the
selected domain type. The adapter checks that two related fields actually select
corresponding SUMMA locations and compatible local axes. A purely numerical
shape check is insufficient when two vectors have different physical meanings.

## Feasibility and application contracts

The central contracts are represented by these existing interfaces:

```fortran
call initialize_summa_calibration_adapter( &
    param_spec, search, param_map, mpar_init, bpar_init, adapter, &
    err=err, message=message)

call is_feasible(adapter%calibration, sampled_values, feasible, &
                 err=err, message=message)
```

`sampled_values` contains untransformed calibration variables, not a flattened
physical field. An ordinary infeasible proposal returns `err = 0` and
`feasible = .false.`. A checker/configuration error returns a nonzero error.
The optimization algorithm owns the retry policy; it does not need to know
which physical relationship failed.

Workers receive accepted scalar vectors. Their adapters rebuild effective fields
from immutable initialized values and write those fields into SUMMA before
parameter-dependent quantities are updated. They do not repeat the coordinator's
full calibration ordered-check pass, although structural/finite safeguards and
SUMMA's model checks remain.

This assumes consistent input initialization on coordinator and workers. If they
read different trial-parameter maps, feasibility checked on one rank may not
describe the model evaluated on another. Data staging is therefore part of the
correctness of the parallel workflow, not just a performance optimization.

## Adding another model

A future model adapter needs to supply parameter metadata and initialized fields,
map generic effective fields into its own storage, and provide model execution
plus observed/simulated-series evaluation. It must also define the meaning of
spatial correspondence and of any local reduction.

The field construction, constraint evaluation, and search proposals can then be
reused without adding model-specific storage into those modules. Do not move
SUMMA-specific NetCDF layouts or native parameter lookups into the generic
feasibility predicate simply to support another model.

Extracting these components into a separate repository/submodule remains a
packaging and integration step. The current small standalone test-model driver
should not be advertised as a complete external-model calibration API.

## Adding a parameterization or objective

A new parameterization needs coordinated updates to parsing, type metadata,
field construction, physical-envelope calculation, application, output
interpretation, and tests. A new ratio reduction must change both construction
and constraint endpoint limits consistently. For example, using local minima in
construction but unrestricted layer extrema in gap calculation would be an
inconsistent definition.

A new target needs a defined observation operator: spatial selection, aggregation,
units, temporal support, missing-value treatment, and metric direction. Naming a
new variable is not enough when its observations require a different vertical or
temporal operator. Update output metadata when introducing a target setting so
results remain interpretable without guessing how the score was produced.

Fortran interfaces keep error outputs last. Optional inputs precede them;
callers that omit such inputs use keyword error arguments. For example:

```fortran
call compute_metric(obs, sim, metric, transform, objective, &
                    boxcox_lambda=0.25_rkind, err=err, message=message)

call compute_metric(obs, sim, metric, transform, objective, &
                    err=err, message=message)
```

