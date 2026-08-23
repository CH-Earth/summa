#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

FILE_MANAGER="$HOME/data/great-slave-lake/Athabasca-river/settings/SUMMA/fileManager_original.txt"

DEV_ROOT="$HOME/models/summa-mpi"
DEV_MPI="$DEV_ROOT/bin/summa.exe"

LOG_DIR="$DEV_ROOT/utils/regression/logs"

echo "File manager:     $FILE_MANAGER"
echo "Dev exe (mpi):    $DEV_MPI"
echo

# ==============================================================================
# Regression runs
# ==============================================================================

mkdir -p "$LOG_DIR"

for np in 1 2 4 8 12 16; do
  echo "Running development MPI: np=$np..."
  mpirun -np "$np" "$DEV_MPI" \
    -m "$FILE_MANAGER" -s "dev_np${np}" \
    > "$LOG_DIR/dev_np${np}.log" 2>&1
done
