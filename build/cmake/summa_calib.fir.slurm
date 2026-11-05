#!/bin/bash
#SBATCH --job-name=summa-calib
#SBATCH --account=rrg-mclark
#SBATCH --nodes=6
#SBATCH --ntasks=1152
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=3:00:00
#SBATCH --output=summa-calib-%j.out
#SBATCH --error=summa-calib-%j.err

set -euo pipefail

# ----------------------------------------------------------------------
# Runtime environment
# ----------------------------------------------------------------------

source "$HOME/models/summa/build/cmake/env_summa.fir.bash"

echo "Job ID:      $SLURM_JOB_ID"
echo "Nodes:       $SLURM_JOB_NODELIST"
echo "CPUs/node:   $SLURM_JOB_CPUS_PER_NODE"
echo "Tasks:       $SLURM_NTASKS"

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
  -n 1152 \
  --ntasks-per-node=192 \
  "$SLURM_TMPDIR/models/summa/bin/summa_sundials_mizuroute_opt.exe" \
  --manifest utils/test/test_calibration/multi_case_example/manifest_century.toml \
  --home "$SLURM_TMPDIR"
