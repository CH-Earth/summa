#!/usr/bin/env bash
#SBATCH --job-name=summa-calib
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# Supply your allocation at submission time: sbatch -A ees250064 ...
# No --mem-per-cpu: Anvil has ~256 GB per standard 128-core node.

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: sbatch -A <account> [--nodes=N] $0 <run_label>" >&2
  exit 2
fi
export RUN_LABEL="$1"

REPO="$HOME/models/summa"
source "$REPO/utils/hpc/anvil/env_summa.anvil.bash"

# The same *path* is created independently on every allocated node.
# Anvil does not automatically set SLURM_TMPDIR.
export SLURM_TMPDIR="${SLURM_TMPDIR:-/tmp/${USER}_${SLURM_JOB_ID}}"

# Shared storage: set SUMMA_DATA_ROOT explicitly if $PROJECT is not
# the allocation project directory containing your Century datasets.
export SUMMA_DATA_ROOT="${SUMMA_DATA_ROOT:-${PROJECT:?PROJECT must be set}/data}"
export SUMMA_OUTPUT_ROOT="${SUMMA_OUTPUT_ROOT:-${SCRATCH:?SCRATCH must be set}/calibration}"

SHARED_EXP="$SUMMA_DATA_ROOT/century/MM/$RUN_LABEL"
if [[ ! -f "$SHARED_EXP/settings/manifest_century.toml" ]]; then
  echo "ERROR: Missing manifest: $SHARED_EXP/settings/manifest_century.toml" >&2
  exit 1
fi

PERSISTENT_DIR="$SUMMA_OUTPUT_ROOT/$RUN_LABEL/$SLURM_JOB_ID"
mkdir -p "$PERSISTENT_DIR"

echo "Job $SLURM_JOB_ID on ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} MPI ranks"
echo "Allocation: ${SLURM_JOB_ACCOUNT:-unknown}; partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "Input root: $SUMMA_DATA_ROOT"
echo "Node-local root: $SLURM_TMPDIR"
echo "Output root: $PERSISTENT_DIR"
module list

# One staging task per node. srun inherits the exported paths/label.
srun --nodes="$SLURM_JOB_NUM_NODES" \
     --ntasks="$SLURM_JOB_NUM_NODES" \
     --ntasks-per-node=1 \
     bash "$REPO/utils/hpc/anvil/setup_summa.century.anvil.bash"

# Identical local manifest path on every node, but separate local filesystems.
LOCAL_MANIFEST="$SLURM_TMPDIR/data/century/MM/$RUN_LABEL/settings/manifest_century.toml"
LOCAL_EXE="$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe"
[[ -f "$LOCAL_MANIFEST" && -x "$LOCAL_EXE" ]] || {
  echo 'ERROR: Staging did not create the local manifest/executable.' >&2
  exit 1
}

# Keep rank groups on each node contiguous (important for cases_per_node).
cd "$SLURM_TMPDIR/models/summa"
srun --distribution=block --kill-on-bad-exit=1 \
     --ntasks="$SLURM_NTASKS" \
     "$LOCAL_EXE" \
     --manifest "$LOCAL_MANIFEST" \
     --home "$SLURM_TMPDIR" \
     --persistent-output "$PERSISTENT_DIR"

# Successful jobs can free local /tmp immediately. A failed job leaves the
# directory for inspection; admins may purge /tmp independently.
if [[ "${KEEP_NODE_SCRATCH:-0}" != "1" ]]; then
  srun --nodes="$SLURM_JOB_NUM_NODES" \
       --ntasks="$SLURM_JOB_NUM_NODES" \
       --ntasks-per-node=1 \
       bash -c '
         expected="/tmp/${USER}_${SLURM_JOB_ID}"
         if [[ "$SLURM_TMPDIR" == "$expected" ]]; then
           rm -rf -- "$SLURM_TMPDIR"
           echo "[$(hostname)] Removed $SLURM_TMPDIR"
         else
           echo "[$(hostname)] Not cleaning unexpected path: $SLURM_TMPDIR" >&2
         fi
       '
fi

echo "Completed calibration: $RUN_LABEL; output: $PERSISTENT_DIR"
