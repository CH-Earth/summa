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

MANIFEST="$HOME/data/century/MM/exp01/settings/manifest_century.toml"

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
echo "Manifest:         $MANIFEST"
echo "--------------------------------------------------"

# ----------------------------------------------------------------------
# Stage SUMMA and data once per node
# ----------------------------------------------------------------------

srun -N "$SLURM_JOB_NUM_NODES" \
     -n "$SLURM_JOB_NUM_NODES" \
     --ntasks-per-node=1 \
     bash -lc 'source $HOME/models/summa/utils/hpc/fir/setup_summa.century.fir.bash'

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
