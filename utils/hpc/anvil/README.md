# Multi-case SUMMA calibration on Anvil

Anvil adaptation of `utils/hpc/fir/`'s Century multi-case calibration scripts.
Check these files into `utils/hpc/anvil/` in the SUMMA repository.

## What changed from Fir

- Full Anvil CPU nodes are **128 cores**, not 192. Use `--ntasks-per-node=128` and
  the `wholenode` partition for production. On `debug`, override this for tests.
- No `--mem-per-cpu=3G`: that exceeds a standard Anvil node's available memory.
- Anvil did not provide `$SLURM_TMPDIR` in our test, so the batch job exports
  `/tmp/${USER}_${SLURM_JOB_ID}` on every compute node. Each node's `/tmp` is separate.
- Stage the executable, domain archive, observations, and **only the selected
  experiment** once on every node. The local manifest and `--home` paths are identical
  across nodes, but resolve to independent node-local filesystems.
- Use Anvil's GCC/OpenMPI/NetCDF/OpenBLAS modules. Build SUNDIALS locally on Anvil
  with a version compatible with the SUMMA code being compiled.
- The script takes an **experiment label** (`exp50`), not a manifest pathname.

## First-time setup (on Anvil)

Check module availability before building (a module shown in documentation may
have changed or require different prerequisites):

```bash
module spider gcc openmpi netcdf-fortran openblas cmake
module spider sundials
module load modtree/cpu
module list
```

The provided environment currently targets the documented GCC 11.2.0 and OpenMPI
4.0.6 stack, with default compatible NetCDF Fortran/OpenBLAS/CMake modules.
If Anvil's `module spider` reports different compatible versions, edit
`env_summa.anvil.bash` **and use the same environment for build and runtime**.

The build script expects `$SUNDIALS_DIR` to be installed (default
`$HOME/local/sundials`). You cannot copy Fir's compiled SUNDIALS binaries.
Use the same SUNDIALS release that works with your Fir SUMMA branch, recompile
it with Anvil's GNU toolchain, and enable the Fortran interface and IDA/IDAS.
For the release-independent installation outline in SUMMA docs, see:
https://summa.readthedocs.io/en/develop_sundials/installation/SUMMA_installation/

Run a lightweight build on a login node or request an interactive debug compute
node for a larger build, then execute:

```bash
cd ~/models/summa/utils/hpc/anvil
bash build.anvil.bash
```

## Shared data layout

Move the inputs from Fir to an Anvil shared filesystem (preferably `$PROJECT`
for retained data). **Do not put large archives in `$HOME`** (25 GB quota).

By default `SUMMA_DATA_ROOT=$PROJECT/data`, with this expected layout:

```text
$SUMMA_DATA_ROOT/
  century/MM/MM-domains.tar                 # must extract domain/<case>/...
  century/MM/exp50/common_inputs/
  century/MM/exp50/settings/manifest_century.toml
  century/MM/exp50/settings/summa_config_template.toml
  camels-spat/observations/obs-daily.tar   # must extract obs-daily/...
```

**Check the archive name before running:** the Fir staging script uses
`MM-domains.tar`, while the supplied Fir README said `MM-exp.tar`. The Anvil
script deliberately follows the *working Fir staging script*. Confirm with
`ls "$SUMMA_DATA_ROOT/century/MM"` and `tar -tf ... | head`.

`$PROJECT` can refer to a different allocation if you have several. Check
`echo "$PROJECT"` or export the exact dataset path in `SUMMA_DATA_ROOT`.
By default finished-basin results go to `$SCRATCH/calibration/<label>/<jobid>`;
Anvil purges scratch based on file age. Set `SUMMA_OUTPUT_ROOT` to a project
path if those results need longer retention.

## Run a test

From the scripts directory (create `logs/` **before** submitting):

```bash
cd ~/models/summa/utils/hpc/anvil
mkdir -p logs
export SUMMA_DATA_ROOT="$PROJECT/data"  # adjust if needed
export SUMMA_OUTPUT_ROOT="$SCRATCH/calibration"

# Short debug run: override ranks per node; update cases_per_node in the manifest
# so the chosen number of ranks is divisible by the number of cases.
sbatch -A ees250064 -p debug --time=00:30:00 \
  --nodes=1 --ntasks-per-node=16 summa_calib.anvil.slurm exp50

# Production (1 full Anvil node = 128 ranks)
sbatch -A ees250064 summa_calib.anvil.slurm exp50

# 4 nodes = 512 ranks; for 4 cases/node: 32 ranks/case, 31 workers
sbatch -A ees250064 --nodes=4 summa_calib.anvil.slurm exp50
```

The manifest's `cases_per_node` must divide the actual ranks per node. Keep
all case groups node-local to avoid stretching one model instance across nodes.
Anvil's current documentation recommends 128-core node-exclusive jobs on the
`wholenode` partition. `srun` without an explicit `--mpi` setting is the
Purdue-documented launch pattern. Before scaling up, test MPI initialization
at low rank counts; if needed inspect `srun --mpi=list` and select the site-
supported PMIx mode.

Monitor `squeue -u "$USER"`, inspect `logs/`, and use `seff <jobid>` after
completion. Local `/tmp` is removed after a successful run unless you submit
with `KEEP_NODE_SCRATCH=1`; on failures it is kept for debugging, but is
still temporary and may be cleared by the cluster.

Anvil references:
- https://docs.rcac.purdue.edu/userguides/anvil/jobs/
- https://docs.rcac.purdue.edu/userguides/anvil/file_management/
- https://docs.rcac.purdue.edu/userguides/anvil/anvil-software/
