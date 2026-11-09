# Optimization algorithms

The SUMMA optimization driver exposes `algorithm = "dds"` and
`algorithm = "nsga2"`. Both operate on the same parameter specification,
physical-field feasibility callback, model evaluation, and MPI trial transport.
Random sampling is used internally for initialization; the generic sampler's
random method is not an additional documented top-level calibration algorithm.

## DDS: one scalar search objective

DDS uses a scalar fitness that is larger for a better fit. When there are
multiple active targets, the scalar is:

```text
F = sum_j [ weight[j] * sign[j] * metric[j] ]

sign = +1 for kge, kgep, nse
sign = -1 for mae, rmse
```

Zero-weight targets are skipped entirely in this arithmetic. Weights must be
finite and nonnegative, with at least one positive weight. With one ordinary
KGE target of weight one, `F` is exactly that KGE value.

There is no automatic normalization by target magnitude, units, variance, or the
sum of weights. Combining RMSE and KGE with equal numerical weights does not
make them equally influential. Define the scalar objective deliberately, or
use separate objectives in NSGA-II when a Pareto interpretation is intended.

### Proposal sequence

DDS fills available workers with feasible random proposals until an elite
individual exists. Thereafter it perturbs the best completed feasible solution.
The probability of selecting each coordinate for perturbation decreases as the
completed-evaluation count approaches the configured budget:

```text
P = 1 - log(max(i, 1)) / log(n_samples)
```

Here `i` is the search progress supplied by the dispatcher. At least one
coordinate is perturbed. The implementation uses a perturbation scale of `0.2`
times that coordinate's transformed search range. Log-transformed parameters
are perturbed in log coordinates and transformed back before feasibility checks.

In this implementation, proposals outside the configured bounds or failing
physical-field constraints are rejected and regenerated. Do not assume another
DDS implementation's boundary reflection rule applies here. Retry limits are
finite.

Workers complete at different times. The coordinator records each result,
updates the elite when appropriate, and supplies the next candidate to an
available worker. It does not wait for a generation boundary. Consequently, a
new proposal can depend on which earlier trials have already returned.

### Budget and practical concurrency

`n_samples` is the number of dispatched model trials per calibration case, not
per rank. Failed model trials use budget; rejected proposals before dispatch do
not. The baseline spin-up and sample-zero export are separate from this budget.

With 191 workers and only 12 samples, most workers cannot evaluate a trial. Even
with a larger budget, many simultaneously outstanding proposals reduce the
number of feedback updates available before subsequent proposals are launched.
More workers can improve throughput without necessarily improving the final
DDS solution at a fixed evaluation budget.

## NSGA-II: separate target objectives

```toml
[calibration]
algorithm = "nsga2"
n_samples = 400

[calibration.nsga2]
population_size = 20
crossover_probability = 0.9
crossover_eta = 20.0
mutation_eta = 20.0
# Omit mutation_probability to use the default based on search dimension.
```

This example uses 20 generations of 20 trials. The sample budget must be an
integer multiple of population size and include at least two generations.
Population size must be at least two.

The first population is randomly generated subject to feasibility. Later
populations use crowded binary tournament parent selection, simulated binary
crossover, and polynomial mutation in transformed search coordinates. Invalid
physical-field candidates are regenerated. Environmental selection retains a
population from parents and offspring using nondomination rank and crowding
distance.

Trials within a generation use the common worker pool, but generation selection
waits for that generation's results. With a population of 20, 191 workers cannot
be kept busy evaluating 191 distinct members of that generation. Group size and
population size should therefore be considered together.

The default population size is 100. The defaults for crossover probability,
crossover distribution index, and mutation distribution index are 0.9, 20, and
20. An omitted mutation probability is derived as one divided by the number of
sampled calibration variables. These settings are not additional model
parameters.

### Directions, weights, and failed evaluations

Efficiency metrics are oriented as maximization objectives and error metrics as
minimization objectives. NSGA-II compares each active objective separately rather
than their DDS weighted sum. The magnitudes of positive weights are ignored.

**Zero-weight targets are rejected for NSGA-II in the current workflow.** There
is no implemented active-objective subset that removes diagnostic columns from
dominance and crowding calculations. Do not set a weight to zero and expect a
quietly ignored Pareto objective.

Failed trials are made unfavorable on every objective. The output includes
population membership and the nondominated flags over the evaluated trial set.
The final nondominated-trial set should not be confused with only the last
retained population. Readers should still exclude failed and incomplete records.

## Independent repeats and random-number initialization

A manifest's `n_runs` creates distinct optimization cases for each basin. The
current case driver initializes the Fortran random-number seed vector using
`42 + run_index - 1`. This is an implementation detail of the present driver,
not a public TOML `seed` option.

A direct `-c` run uses the default run index of one. Repeating that command, or
changing only an output suffix, is not a controlled way to request a different
seed. Use manifest-driven repeats when independent run indices are required.

Matching run indices across experiments is useful organizationally, but does
not guarantee matched parameter proposals. Search dimension, feasible region,
completed-trial order, compiler/runtime random-number implementation, and
numerical behavior can differ.

For asynchronous DDS, even identical initial seeds and model inputs do not
promise identical trajectories across worker counts or machines. To test model
numerical consistency, compare **identical physical parameter fields and initial
states**, not the final winners of independent asynchronous searches.

## Restart, replay, and continuation are different

A shared spin-up restart initializes every trial's model state. A selected
physical parameter export supports a standalone simulation. Neither contains
the complete state of the optimization algorithm.

The current driver does not read a prior calibration NetCDF to restore DDS elite
state, pending tasks, random-number state, or NSGA-II population/generation state.
Increasing `n_samples` and launching again starts another search, not a supported
continuation of the old one. Changing the budget also changes the DDS progress
schedule.

Optimization checkpoints, population seeding from a coarse-resolution
calibration, and recovery of in-flight MPI work require explicit additional
implementation. Treat them as future extensions, not guarantees provided by
retaining model restart files.

**Implementation entry points:** `summa_parameter_sampling.f90`,
`parameter_search.f90`, `nsga2.f90`, and `model_calibration.f90` in
`build/source/calibration/`; `scalarize_objectives` and `oriented_objectives` in
`build/source/driver/summa_simulation.f90`; and `run_case` in
`build/source/driver/summa_calibration_case.f90`.
