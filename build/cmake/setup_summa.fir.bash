#!/bin/bash

# SUMMA runtime environment on Fir

# ----------------------------------------------------------------------
# Load modules
# ----------------------------------------------------------------------

module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load openblas/0.3.24
module load netcdf-mpi/4.9.2
module load netcdf-fortran-mpi/4.6.1

# SUNDIALS installation
export SUNDIALS_DIR=$HOME/local/sundials

# ----------------------------------------------------------------------
# Limit each MPI model instance to one compute thread
# ----------------------------------------------------------------------

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

# ----------------------------------------------------------------------
# Stage SUMMA to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/models/summa"

# Executables
cp -a "$HOME/models/summa/bin" \
      "$SLURM_TMPDIR/models/summa/"

# Multi-case configuration and test utilities
mkdir -p "$SLURM_TMPDIR/models/summa/utils/test/test_calibration"

cp -a "$HOME/models/summa/utils/test/test_calibration" \
      "$SLURM_TMPDIR/models/summa/utils/test"

# ----------------------------------------------------------------------
# Stage model inputs and observations to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/century/test"
mkdir -p "$SLURM_TMPDIR/data/camels-spat/observations"

# Copy archives as single files to minimize Lustre metadata I/O

echo "[$(hostname)] Copying Century archive..."
cp "$HOME/data/century/test/exp01.tar" \
   "$SLURM_TMPDIR/data/century/test/"
echo "[$(hostname)] Century archive copied"

echo "[$(hostname)] Copying observations archive..."
cp "$HOME/data/camels-spat/observations/obs-daily.tar" \
   "$SLURM_TMPDIR/data/camels-spat/observations/"
echo "[$(hostname)] Observations archive copied"

# Extract archives from node-local storage

echo "[$(hostname)] Extracting Century archive..."
tar --warning=no-unknown-keyword \
    -xf "$SLURM_TMPDIR/data/century/test/exp01.tar" \
    -C "$SLURM_TMPDIR/data/century/test/"
echo "[$(hostname)] Century archive extracted"

echo "[$(hostname)] Extracting observations archive..."
tar --warning=no-unknown-keyword \
    -xf "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar" \
    -C "$SLURM_TMPDIR/data/camels-spat/observations/"
echo "[$(hostname)] Observations archive extracted"

# Create case-specific runtime output directories
for case_dir in "$SLURM_TMPDIR"/data/century/test/exp01/domain/*; do
  mkdir -p "$case_dir/work"
done

# Remove local archive copies after extraction
rm "$SLURM_TMPDIR/data/century/test/exp01.tar"
rm "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar"

# ----------------------------------------------------------------------
# Move to node-local working directory
# ----------------------------------------------------------------------

cd "$SLURM_TMPDIR/models/summa"

echo "SUMMA runtime environment ready"
echo "Host: $(hostname)"
echo "Local work directory: $SLURM_TMPDIR/summa"
echo "Local data directory: $SLURM_TMPDIR/data"
echo "Available CPUs: $SLURM_CPUS_ON_NODE"
