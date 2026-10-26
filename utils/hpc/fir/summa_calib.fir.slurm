#!/bin/bash
#SBATCH --job-name=summa-calib
#SBATCH --account=rrg-mclark
#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=summa-calib-%j.out
#SBATCH --error=summa-calib-%j.err

# Multi-case calibration of the Century basins on Fir.
#   sbatch [--nodes=N --ntasks=192*N] summa_calib.fir.slurm [manifest]
# The manifest's cases_per_node sets how many basins share a node. PR 649 scaling runs:
#   --nodes=1,2,4,8 with cases_per_node=4, and --nodes=2 with cases_per_node=1,2,4,8

set -euo pipefail

MANIFEST=${1:-"$HOME/data/century/MM/exp01/settings/manifest_century.toml"}

source "$HOME/models/summa/utils/hpc/fir/env_summa.fir.bash"

echo "Job ID:      $SLURM_JOB_ID"
echo "Nodes:       $SLURM_JOB_NODELIST"
echo "Tasks:       $SLURM_NTASKS"
echo "Manifest:    $MANIFEST"

# stage SUMMA and data once per node
srun -N "$SLURM_JOB_NUM_NODES" \
     -n "$SLURM_JOB_NUM_NODES" \
     --ntasks-per-node=1 \
     bash -lc 'source $HOME/models/summa/utils/hpc/fir/setup_summa.century.fir.bash'

cd "$SLURM_TMPDIR/models/summa"

srun --mpi=pmix \
  -n "$SLURM_NTASKS" \
  --ntasks-per-node=192 \
  "$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe" \
  --manifest "$MANIFEST" \
  --home "$SLURM_TMPDIR" \
  --persistent-output "$SCRATCH/calibration/$SLURM_JOB_ID"
