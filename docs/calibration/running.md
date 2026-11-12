# Configuration and execution

Start with a working simulation and retain its physics, attributes, routing
configuration, forcing organization, and observation conventions. Calibration
adds a search specification; it does not create the underlying domain inputs.
The examples below are structural examples, not self-contained test datasets.

## Build the optimization executable

A CMake build with `USE_MPI=ON` produces the optimization executable. Its name
reflects enabled components. For example, a build with SUNDIALS and mizuRoute
uses `summa_sundials_mizuroute_opt.exe`; a different build may omit one or both
suffixes.

From the repository root, an illustrative configuration is:

```bash
cmake -S build -B build/cmake_build \
  -DUSE_MPI=ON \
  -DUSE_SUNDIALS=ON \
  -DUSE_MIZUROUTE=ON \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build/cmake_build --parallel 8
```

First supply the compiler, MPI, NetCDF-Fortran, SUNDIALS, and other dependencies
required by your installation. Use the repository's platform-specific build
scripts where appropriate. See [Installation](../installation/SUMMA_installation.md)
and [mizuRoute build instructions](../mizuroute/build.md).

The finite-value checks in calibration and metric evaluation depend on normal
IEEE behavior. Avoid compiler settings that assume NaNs/infinities cannot exist,
such as finite-math-only assumptions. Rebuild dependent Fortran modules when
interfaces change, and reconfigure CMake when source lists or build-time version
information need refreshing.

## A direct single-case configuration

[dds_case.toml](examples/dds_case.toml) supplies a complete structural example
using an uncoupled domain, three physical calibration variables, and KGE plus
two zero-weight diagnostics. It expects files such as:

```text
ROOT/
  data/example/
    domain/CAN_01AD003/
      summa_inputs/       attributes, forcing list, trial parameters
      summa_state/        coldState.nc and generated restart state
      summa_forcing/      forcing NetCDF files
    exp01/
      common_inputs/      decisions, output control, parameter tables
      work/               writable experiment-specific model output
    observations/         observed time series
```

Copy the example and replace paths, filenames, parameter bounds, dates, and
observation choices with those for your working domain. Its small sample count
is for an initial execution check, not a scientifically meaningful calibration.

```bash
mpirun -np 12 ./bin/summa_sundials_mizuroute_opt.exe \
  -c /absolute/path/to/dds_case.toml \
  --home /absolute/path/to/ROOT \
  --persistent-output /absolute/path/to/results/exp01/job-test
```

There is one search and eleven workers, not twelve independent optimization
runs. The direct example uses a literal experiment directory in its paths.

For coupled routing, retain the working `[mizuRoute]`, `[hydrofabric]`, and any
required `[remapping]` sections, set `use_mizuroute = true`, and select the
appropriate observed reach in each target. Copying the calibration section
alone does not replace the routing setup. See [mizuRoute execution](../mizuroute/running.md).

## Path resolution

`--home` supplies the root substituted for `{home}`. It does not automatically
copy files or rewrite absolute paths embedded inside forcing lists, routing
namelists, or symbolic links.

In the illustrated layout:

```toml
[simulation]
case_name = "CAN_01AD003"
basin_dir = "data/example/domain/{case_name}"

[summa_files]
settings_path = "{home}/"
state_path = "{home}/{basin_dir}/summa_state/"
forcing_path = "{home}/{basin_dir}/summa_forcing/"
output_path = "{home}/data/example/exp01/work/{case_name}/"
init_condition = "coldState.nc"
attributes = "{basin_dir}/summa_inputs/attributes.nc"
```

`attributes`, the forcing-list filename, trial parameters, decisions, and table
filenames resolve relative to `settings_path`. The initial-condition filename
resolves relative to `state_path`. The state, forcing, and output paths themselves
are expanded filesystem paths. Preserve trailing directory separators because
several model interfaces concatenate a directory and filename directly.

A case name must match the directory naming convention. A directory `07DD001`
and a case label `CAN_07DD001` do not become identical automatically. Either
use the matching case name or set `basin_dir` explicitly. The observation filename
may legitimately retain a different prefix if it is specified explicitly.

## Manifest, template, and experiment tag

[manifest.toml](examples/manifest.toml) and
[dds_template.toml](examples/dds_template.toml) demonstrate two basins, two
independent runs per basin, and two case groups per node:

```toml
[multi_case]
experiment = "exp01"
template_file = "dds_template.toml"
case_names = ["CAN_01AD003", "CAN_01AM001"]
n_runs = 2
cases_per_node = 2
```

`template_path` is optional and defaults to the manifest's directory. Set it
explicitly when the template lives elsewhere. The template uses:

```toml
[simulation]
case_name = "{case_name}"
basin_dir = "data/example/domain/{case_name}"

[summa_files]
settings_path = "{home}/"
decisions = "data/example/{experiment}/common_inputs/modelDecisions.txt"
output_path = "{home}/data/example/{experiment}/work/{case_name}/"
```

These are excerpts from complete files, not additional duplicate TOML sections
to append to an existing configuration.

Launch the manifest rather than its unresolved template:

```bash
mpirun -np 12 ./bin/summa_sundials_mizuroute_opt.exe \
  --manifest /absolute/path/to/manifest.toml \
  --home /absolute/path/to/ROOT \
  --persistent-output /absolute/path/to/results/exp01/job-test
```

On one 12-rank node with two groups, this provides two concurrent cases with
five trial workers each. Four total cases are processed as groups become free.

**A manifest can contain one basin.** Use this to obtain an experiment tag and
independent repeats for a distributed basin without maintaining a separate
single-basin scheduling mechanism.

In the checked implementation, the experiment tag is parsed from
`[multi_case]`; `[simulation].experiment` is not a supported parser entry.
For direct `-c` configurations, use literal experiment paths, or use a one-basin
manifest when `{experiment}` is needed. Do not infer that a directory named
`exp55` automatically defines the placeholder.

## Spin-up and scoring periods

The optimization driver performs a baseline spin-up for the year preceding
`simulation.start_time`. It derives that interval by moving the start year back
one year. For the examples, forcing must therefore extend back to October 2014,
not just the 2015 trial start.

Each rank obtains baseline initialized parameter fields. The group leader writes
the common spin-up state, and each trial begins from that state using its proposed
parameters. The spin-up is not rerun to equilibrium for every candidate.

There are three distinct periods:

| Period | Purpose |
| --- | --- |
| Baseline spin-up before simulation start | Produce the common initial state using the baseline configuration |
| Trial simulation before scoring starts | Permit candidate-specific adjustment without scoring those observations |
| Scoring period | Evaluate aligned observations selected by the calibration dates |

A later `calibration.start_date` does not shorten the simulated period; it changes
which observations are scored. For slowly adjusting aquifer storage, test the
sensitivity to the shared initial state and the unscored trial interval before
interpreting parameter differences scientifically.

The date-window semantics for period-ending observations are described in
[Objectives](objectives.md). The simulation end and observation filter end are
not interchangeable notions of an inclusive final day.

## Interactive Slurm test

Request resources using your site's policy, then stage data on the assigned
compute node. A staging script invoked as a separate `srun` step changes that
step's environment, not the login shell's environment. Resolve node-local
variables inside a compute-node shell:

```bash
srun --nodes=1 --ntasks=1 bash -lc '
  source "$HOME/models/summa/utils/hpc/fir/setup_summa.distributed.fir.bash"
'

srun --mpi=pmix --nodes=1 --ntasks=12 bash -lc '
  source "$HOME/models/summa/utils/hpc/fir/env_summa.fir.bash"
  exec "$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe" \
    -c "$SLURM_TMPDIR/data/distributed/FH/exp55/settings/summa_config_template.toml" \
    --home "$SLURM_TMPDIR" \
    --persistent-output "$SCRATCH/calibration/FH-exp55/$SLURM_JOB_ID"
'
```

This illustrates the Fir distributed-domain layout used during development. It
requires that staging script and the named configuration to exist in your
checkout/data installation. Other sites need their own environment and launcher.
The allocation must already contain at least the requested tasks and memory.

For a smoke test, change only the staged configuration to a small positive DDS
budget that allows returning workers to be redispatched. Verify paths, parameter
application, observations, outputs, and persistence. Do not use that tiny budget
to rank scientific parameterizations. Restaging can overwrite the edited
node-local configuration, so make the sequence explicit.

## Multi-node production staging

For manifest jobs, run staging once on every node before launching all ranks:

```bash
# Inside a batch allocation, after loading the site's environment.
srun --nodes="$SLURM_JOB_NUM_NODES" \
     --ntasks="$SLURM_JOB_NUM_NODES" \
     --ntasks-per-node=1 bash -lc '
  source "$HOME/models/summa/utils/hpc/fir/setup_summa.century.fir.bash"
'

srun --mpi=pmix --ntasks="$SLURM_NTASKS" bash -lc '
  source "$HOME/models/summa/utils/hpc/fir/env_summa.fir.bash"
  exec "$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe" \
    --manifest "$SLURM_TMPDIR/data/century/MM/exp55/settings/manifest_century.toml" \
    --home "$SLURM_TMPDIR" \
    --persistent-output "$SCRATCH/calibration/exp55/$SLURM_JOB_ID"
'
```

The task placement in the allocation should be uniform, for example equal
ranks per node with one CPU per rank for an MPI-only run. The job must request
memory appropriate to the measured workload. Neither the executable nor
`--home` requests additional scheduler resources.

Keep forcing and observations separate from small experiment-specific settings.
A reusable archive can contain `domain/`, while experiment `common_inputs/` and
`settings/` are refreshed independently. Exclude runtime work directories from
the domain archive, check external symlinks, and create writable work locations
after extraction. Do not rerun a script that deletes/refills shared runtime
folders while another step in that allocation is using them.

