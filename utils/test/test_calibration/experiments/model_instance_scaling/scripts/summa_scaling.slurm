#!/bin/bash
#SBATCH --job-name=summa-scaling
#SBATCH --account=rrg-mclark
#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# ----------------------------------------------------------------------
# Experiment configuration
# ----------------------------------------------------------------------

if [ "$#" -ne 1 ]; then
  echo "Usage: sbatch [options] $0 <cases_per_node>"
  echo "Example: sbatch --nodes=2 --ntasks=384 $0 4"
  exit 1
fi

CASES_PER_NODE="$1"

case "$CASES_PER_NODE" in
  1|2|4|8)
    ;;
  *)
    echo "ERROR: cases_per_node must be one of: 1, 2, 4, 8"
    exit 1
    ;;
esac

RANKS_PER_NODE=192
RANKS_PER_CASE=$((RANKS_PER_NODE / CASES_PER_NODE))
WORKERS_PER_CASE=$((RANKS_PER_CASE - 1))

MANIFEST="utils/test/test_calibration/experiments/model_instance_scaling/configs/manifest_c${CASES_PER_NODE}.toml"

# ----------------------------------------------------------------------
# Runtime environment
# ----------------------------------------------------------------------

source "$HOME/models/summa/build/cmake/env_summa.fir.bash"

echo "--------------------------------------------------"
echo "SUMMA model-instance scaling experiment"
echo "--------------------------------------------------"
echo "Job ID:           $SLURM_JOB_ID"
echo "Job name:         $SLURM_JOB_NAME"
echo "Nodes:            $SLURM_JOB_NUM_NODES"
echo "Node list:        $SLURM_JOB_NODELIST"
echo "Tasks:            $SLURM_NTASKS"
echo "Cases/node:       $CASES_PER_NODE"
echo "Ranks/case:       $RANKS_PER_CASE"
echo "Workers/case:     $WORKERS_PER_CASE"
echo "Manifest:         $MANIFEST"
echo "--------------------------------------------------"

# ----------------------------------------------------------------------
# Stage SUMMA and data once per node
# ----------------------------------------------------------------------

srun -N "$SLURM_JOB_NUM_NODES" \
     -n "$SLURM_JOB_NUM_NODES" \
     --ntasks-per-node=1 \
     bash -lc 'source $HOME/models/summa/build/cmake/setup_summa.fir.bash'

# ----------------------------------------------------------------------
# Run calibration
# ----------------------------------------------------------------------

cd "$SLURM_TMPDIR/models/summa"

srun --mpi=pmix \
  -n "$SLURM_NTASKS" \
  --ntasks-per-node=192 \
  "$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe" \
  --manifest "$MANIFEST" \
  --home "$SLURM_TMPDIR" \
  --persistent-output "$SCRATCH/calibration/$SLURM_JOB_ID"
