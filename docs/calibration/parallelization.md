# Parallelization strategy

The calibration workflow uses two scheduling levels: independent calibration
cases, and independent parameter trials within a case. Model-domain
decomposition is a separate concept and is not currently combined with those
levels inside the optimization driver.

## Three kinds of parallel work

| Level | Work being divided | Current calibration behavior |
| --- | --- | --- |
| Cases | Different basins or independent optimization repeats | Manifest cases are dynamically assigned to node-local case groups |
| Model instances | Different parameter trials for one case | One coordinator dispatches trials to independent MPI workers |
| Model domain | GRUs/HRUs within one parameter trial | Not decomposed across trial workers; each trial uses `MPI_COMM_SELF` |

A 69-HRU Athabasca configuration remains a 69-HRU simulation on **each evaluation
worker**. It is not divided into 69 tasks merely because the job has many MPI
ranks. Each worker has its own model state, forcing access, routing instance,
and parameter application for the candidate it is evaluating.

## Direct single-case execution: `-c`

Without a manifest, the instance communicator is the world communicator. All
launched ranks participate in **one** calibration search:

```text
MPI_COMM_WORLD = one calibration instance group

rank 0: coordinator, search state, feasibility, output
  |
  +-- rank 1: complete model evaluation for sample a
  +-- rank 2: complete model evaluation for sample b
  +-- ...
  +-- rank R-1: complete model evaluation for another sample
```

On one node, `R = 12` gives eleven workers. `R = 192` gives 191 workers. The
coordinator is part of the requested rank count; it is not an additional rank
outside the allocation.

There is no node-local `cases_per_node` partition in this direct path. If the
world communicator spans several nodes, it is still one search pool, not one
independent calibration on each node. That arrangement requires shared access
to its runtime files, especially the common spin-up restart; see below.

## Manifest execution: `--manifest`

The manifest identifies a template, a list of case names, `n_runs`, and
`cases_per_node`. The work list contains:

```text
number of cases = number of case names * n_runs
```

For two basins and two repeats, the logical work items are basin A/run1,
basin A/run2, basin B/run1, and basin B/run2. Groups claim them as capacity
becomes available; a basin is not permanently assigned to a particular node.

The parallel setup uses a shared-memory split of the world communicator to
identify ranks on each node. It then divides each node into case groups. A
node-leader communicator establishes node numbering; the resulting case-group
communicators perform calibration work.

For a uniform allocation with `R` ranks per node and `C` case groups per node:

```text
ranks per case group   = R / C
workers per case group = R / C - 1
maximum active cases   = number of nodes * C
```

Every group needs at least a coordinator and one worker. The current setup
reduces an excessive requested group count and reduces a non-dividing count
until it divides the node rank count. Request a valid layout deliberately and
inspect the printed rank/group mapping rather than relying on this adjustment.
Use the same ranks per node throughout a manifest job; group numbering and
reported totals assume the homogeneous layout used by these examples.

### Example: one 192-rank node

| Cases per node | Ranks per case | Workers per case | Total evaluation workers |
| --- | ---: | ---: | ---: |
| 1 | 192 | 191 | 191 |
| 2 | 96 | 95 | 190 |
| 4 | 48 | 47 | 188 |
| 8 | 24 | 23 | 184 |

With four groups, a node looks like:

```text
node
  case group 0: coordinator + 47 workers -> one basin/run
  case group 1: coordinator + 47 workers -> another basin/run
  case group 2: coordinator + 47 workers -> another basin/run
  case group 3: coordinator + 47 workers -> another basin/run
```

Each group drains its current case and then claims another. Cores are not
reassigned between groups in the middle of a case, and idle workers are not
borrowed by a different case group.

### Dynamic assignment across nodes

A global case counter is exposed through an MPI window. Each case-group leader
claims the next case, broadcasts the index to its own workers, and invokes the
model-specific case callback. When that case finishes, the group requests the
next one. A zero assignment signals that no more cases remain.

This avoids a fixed list of basins per node when basin runtimes differ. It does
not eliminate the final job tail once no cases remain to replace finishing work.
There is no separate global scheduler rank that needs to be added to the rank
counts in the table.

Manifest groups are node-local. Adding nodes adds more independent groups; it
does not make one group or one model evaluation span those nodes.

## Trial dispatch within a case group

The coordinator retains the optimization state, checks candidate feasibility,
and sends a sample identifier and the scalar calibration vector to a worker.
It does not send a complete HRU-by-layer parameter map for every proposal.
Each worker's model-specific adapter constructs that map from its retained
initialized fields and the received coefficients.

A worker initializes SUMMA from the case's shared spin-up state, applies the
candidate, runs the complete domain, and returns the objective vector plus the
trial-failure status. The coordinator writes the trial record and dispatches
more work.

DDS uses arriving results immediately. NSGA-II shares the same evaluation
transport but waits at generation boundaries before selecting the next
population. These synchronization rules are why the same worker allocation
can have different efficiency for the two algorithms.

## Why use several cases per node?

More cases per node allocate fewer workers to any one search and add coordinator
overhead. However, they can reduce the impact of a small number of slow trials
on overall throughput across a large basin collection. They also allow a search
to receive more sequential feedback relative to the number of outstanding
proposals at a fixed budget.

There is no universally optimal case-group count. Compare completed evaluations
per unit wall time, total case completion time, idle tails, memory, and achieved
calibration performance. Do not infer that doubling workers halves runtime or
improves the best score at a fixed budget.

For NSGA-II, a group's useful concurrent workers are also limited by the number
of candidates in a generation. For DDS, a short smoke test should have enough
trials to exercise returning workers and redispatch; one initial wave is not a
scaling benchmark.

## Data locality and the shared restart

Every rank in a case performs baseline initialization and spin-up to obtain the
initialized parameter structures. The case leader writes the shared spin-up
restart; its filename is broadcast to the group. The state file's contents are
not broadcast as an MPI memory object.

For a node-local manifest group, all group members can read that node's local
restart. Stage the domain once per node, not once per worker. On different nodes,
identically spelled `$SLURM_TMPDIR` paths do not imply a shared filesystem.

For a direct `-c` pool spanning nodes, all workers must be able to read the
restart written by rank zero. Shared storage or an explicitly implemented
replication workflow is necessary. Simply staging independent input copies on
every node does not distribute the subsequently generated restart.

Use distinct output locations for different experiments and jobs. Within one
calibration, run/sample/rank naming separates work; it does not protect two
separate jobs intentionally writing the same experiment, case, and run into the
same output directory.

## Memory and threading

Each worker holds a complete model instance. Memory is therefore driven by
per-instance storage multiplied by active workers, plus coordinator and runtime
overheads. More allocated cores do not reduce the domain size held by one
calibration worker.

For a pure MPI baseline, set library/OpenMP thread counts to one unless a tested
hybrid arrangement is intended:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
```

For large or long distributed runs, consider the runtime model decisions:

```text
read_force  readPerStep
write_buff  writePerStep
```

These limit full-series forcing/output buffering. They do not eliminate retained
calibration series, routing arrays, or other allocations. `write_timeseries =
false` disables normal model time-series writing but is not equivalent to choosing
per-step buffering. Measure memory after spin-up and during complete trials,
not only immediately after launching.

Per-step reading can increase I/O calls. Node-local staging can help separate
shared-filesystem traffic from model execution, but performance remains an
empirical property of the domain, forcing-file organization, and worker count.

## Current boundary and future extensions

The `_mpi.exe` simulation executable partitions model GRUs; the `_opt.exe`
calibration executable distributes parameter trials. Building with both MPI
and mizuRoute enables parallel **independent coupled instances**, not MPI routing
or domain decomposition within each trial.

Combining multi-rank domain evaluation with multiple concurrent instances would
require another communicator hierarchy, collective evaluation/parameter
application within each trial, and compatible routing and file handling. The
current modular structure helps organize that extension but does not implement
it automatically.

**Implementation entry points:** `summa_driver_opt.f90` and
`summa_calibration_case.f90` in `build/source/driver/`;
`calibration_parallel.f90`, `dynamic_case_dispatch.f90`, and
`summa_parameter_sampling.f90` in `build/source/calibration/`.
