#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

FILE_MANAGER="$HOME/data/great-slave-lake/Athabasca-river/settings/SUMMA/fileManager.txt"

DEV_ROOT="$HOME/models/summa"
STABLE_ROOT="$HOME/models/stable/summa"

DEV_SERIAL="$DEV_ROOT/bin/summa_sundials.exe"
DEV_MPI="$DEV_ROOT/bin/summa_sundials_mpi.exe"
STABLE_SERIAL="$STABLE_ROOT/bin/summa_sundials.exe"

LOG_DIR="$DEV_ROOT/utils/regression/logs"

echo "File manager:     $FILE_MANAGER"
echo "Stable exe:       $STABLE_SERIAL"
echo "Dev exe (serial): $DEV_SERIAL"
echo "Dev exe (mpi):    $DEV_MPI"
echo

# ==============================================================================
# Regression runs
# ==============================================================================

mkdir -p "$LOG_DIR"

echo "Running stable serial..."
"$STABLE_SERIAL" -m "$FILE_MANAGER" -s stable_serial > "$LOG_DIR/stable_serial.log" 2>&1

echo "Running development serial..."
"$DEV_SERIAL"    -m "$FILE_MANAGER" -s dev_serial    > "$LOG_DIR/dev_serial.log"    2>&1

for np in 1 2 4 8 12 16; do
  echo "Running development MPI: np=$np..."
  mpirun -np "$np" "$DEV_MPI" \
    -m "$FILE_MANAGER" -s "dev_np${np}" \
    > "$LOG_DIR/dev_np${np}.log" 2>&1
done
