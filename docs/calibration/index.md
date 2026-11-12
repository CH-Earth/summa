# Parameter calibration

SUMMA's calibration driver repeatedly initializes and runs a model, compares
simulated series with observations, and searches for useful parameter values.
It supports both collections of independent catchments and individual spatially
distributed model domains, with optional in-process mizuRoute coupling.

The implementation separates three concerns:

1. **Search and feasibility:** propose scalar calibration variables, construct
   their effective physical fields, and check configured constraints.
2. **Model evaluation:** initialize SUMMA, apply those fields, run the model,
   collect simulated series, and calculate target metrics.
3. **Parallel execution:** schedule independent cases and parameter trials using
   MPI, and record completed evaluations.

This separation is important. An optimization algorithm does not need to know
how SUMMA stores a layered soil parameter, and a parameter adjustment does not
need to know whether DDS or NSGA-II proposed it.

## Guide

| Page | Contents |
| --- | --- |
| [Parameter adjustments](parameters.md) | Physical values, spatial multipliers, tied ratios, chains, and search transformations |
| [Constraints and feasibility](constraints.md) | Physical ordering, minimum gaps, spatial compatibility, and rejection sampling |
| [Objectives and observations](objectives.md) | Targets, spatial selection, time alignment, metrics, transformations, and diagnostic targets |
| [Optimization algorithms](algorithms.md) | DDS, NSGA-II, budgets, seeds, reproducibility, and restart limitations |
| [Parallelization](parallelization.md) | Cases, instance groups, MPI workers, node-local scheduling, and domain-decomposition limits |
| [Configuration and execution](running.md) | Direct and manifest runs, examples, spin-up, staging, memory, and Slurm |
| [Outputs and diagnostics](outputs.md) | NetCDF contents, selecting parameters, physical exports, provenance, and persistence |
| [Implementation and extension](implementation.md) | Module responsibilities, execution sequence, adapter contracts, and testing |
| [Troubleshooting](troubleshooting.md) | Preflight checks and common failure modes |

## Terms used throughout this guide

A **simulation case** is one independently configured and executed SUMMA
simulation and, when calibration is enabled, its associated optimization run.
A simulation case can contain multiple grouped response units (GRUs),
hydrologic response units (HRUs), and sub-HRU domains. A distributed river
basin is treated as one simulation case, regardless of its spatial complexity.
Multiple simulation cases can represent the same basin, for example when
performing independent calibration repeats with different random seeds. In
a manifest, two basins with two repeats each produce four simulation cases.

A **trial** is one model evaluation of the simulation case using a proposed parameter vector. One trial
can return several objective values. A **target** specifies the observations,
simulated variable, spatial selection, metric, and transformation for one such
value.

A **calibration variable** is a scalar searched by the optimizer. A **physical
parameter field** is the collection of parameter values actually applied to the
model. They are identical concepts only for a uniform physical-value adjustment.
A single multiplier can control a field containing thousands of elements.

A **case group** is an MPI communicator containing one coordinating rank and
its evaluation workers. A **model instance** is a worker's initialized SUMMA
simulation, not another independently calibrated basin.

## What is implemented

The current workflow supports:

- Scalar physical-value adjustments, scaling of initialized spatial fields, and
  ratios based on effective reference fields, including chained dependencies.
- A shared physical-field feasibility predicate used by the search algorithms.
- DDS with a weighted scalar objective, and NSGA-II with separate objectives.
- Additional zero-weight diagnostic targets in DDS.
- Concurrent trials, and manifest-driven scheduling of independent basin/run
  cases across node-local MPI groups.
- Calibration NetCDF output, selected physical parameter exports, and a small
  provenance record.

**Each trial currently runs its complete model domain on one MPI rank.**
The optimization driver creates a single-rank domain communicator for that
trial. Giving a calibration n ranks provides n-1 concurrent evaluation
workers; it does not partition one trial's HRUs across those workers. The
separate domain-parallel simulation executable is a different execution path.
See [Parallelization](parallelization.md).

## What is not implied by these capabilities

The current driver does not combine MPI domain decomposition inside a trial
with MPI parallelism across trials. Manifest case groups do not span nodes.
A direct single-case worker pool can span nodes, but its shared restart and input
paths must then be accessible to every participating rank.

A model-state restart is not an optimizer checkpoint. Calibration output records
trials, but the driver does not restore the DDS or NSGA-II search state from that
file after interruption. Coarse-to-fine population seeding and resumable
optimization are separate extensions, not implicit behavior of parameter export.

The generic calibration components remain part of the SUMMA repository. Their
separation prepares them for reuse by other models; it is not yet a separately
packaged, model-independent calibration service.

## Start from a working simulation

First establish a working SUMMA configuration, including the observation time
convention and any routing inputs. Then add calibration settings, parameters,
and targets. Run a small calibration before increasing the budget or concurrency.

For a build with MPI, SUNDIALS, and mizuRoute enabled, a direct run is:

```bash
mpirun -np 12 ./bin/summa_sundials_mizuroute_opt.exe \
  -c /path/to/case.toml \
  --home /path/to/staged/root
```

This is one calibration with one coordinator and eleven trial workers. Coupled
routing is activated separately by `simulation.use_mizuroute = true`; an executable
built with routing support can also run an uncoupled configuration.

Use a manifest for independent repeats, multiple basins, or multiple case groups
per node. A manifest may contain just one basin. The [execution guide](running.md)
provides complete structural examples and explains which paths must be changed
for a real dataset.

