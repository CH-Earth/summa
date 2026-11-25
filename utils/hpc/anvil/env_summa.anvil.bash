#!/usr/bin/env bash
# Source this file to select the Anvil CPU toolchain for SUMMA.
# Versions should be verified using `module spider` on Anvil before first use.

if ! type module >/dev/null 2>&1; then
  echo 'ERROR: Lmod `module` command is unavailable. Start an Anvil Bash login shell.' >&2
  return 1 2>/dev/null || exit 1
fi

module purge
module load modtree/cpu
module load gcc/11.2.0
module load openmpi/4.0.6
module load netcdf-fortran
module load openblas

# Use locally installed CMake (SUMMA requires >= 3.21)
export PATH="$HOME/local/cmake-3.28.4-linux-x86_64/bin:$PATH"

# Build SUNDIALS on Anvil with the same tested release as on Fir.
export SUNDIALS_DIR="${SUNDIALS_DIR:-$HOME/local/sundials}"

# Prevent one model instance from spawning multiple threaded BLAS/OpenMP workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
