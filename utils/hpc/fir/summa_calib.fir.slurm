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
# Command-line arguments
# ----------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <experiment_directory> <setup_script>"
  exit 1
fi

EXPERIMENT_DIR="$(realpath "$1")"
SETUP_SCRIPT="$(realpath "$2")"

if [ ! -d "$EXPERIMENT_DIR" ]; then
  echo "ERROR: Experiment directory does not exist: $EXPERIMENT_DIR" >&2
  exit 1
fi

if [ ! -f "$SETUP_SCRIPT" ]; then
  echo "ERROR: Setup script does not exist: $SETUP_SCRIPT" >&2
  exit 1
fi

# ----------------------------------------------------------------------
# Experiment configuration
# ----------------------------------------------------------------------

# Derive the run label from the experiment directory name.
RUN_LABEL="$(basename "$EXPERIMENT_DIR")"

MANIFEST="$EXPERIMENT_DIR/settings/manifest.toml"

if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: Manifest does not exist: $MANIFEST" >&2
  exit 1
fi

echo "Experiment directory: $EXPERIMENT_DIR"
echo "Run label:            $RUN_LABEL"
echo "Manifest:             $MANIFEST"
echo "Setup script:         $SETUP_SCRIPT"
echo

# ----------------------------------------------------------------------
# Runtime environment
# ----------------------------------------------------------------------

source "$HOME/models/summa/utils/hpc/fir/env_summa.fir.bash"

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
     bash -lc "source '$SETUP_SCRIPT'"

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
  --persistent-output "$SCRATCH/calibration/${RUN_LABEL}/${SLURM_JOB_ID}"
