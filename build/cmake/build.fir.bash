#!/usr/bin/env bash
set -e

# --------------------------------------------------------------------------------------
# SUMMA build environment for fir
# --------------------------------------------------------------------------------------

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load netcdf-mpi/4.9.2
module load netcdf-fortran-mpi/4.6.1
module load openblas/0.3.24

# SUNDIALS installation
export SUNDIALS_DIR="$HOME/local/sundials"

# Avoid inode exhaustion in /tmp
mkdir -p "$SCRATCH/tmp"
export TMPDIR="$SCRATCH/tmp"

# --------------------------------------------------------------------------------------
# Configure and build SUMMA
# --------------------------------------------------------------------------------------

cd "$HOME/models/summa/build/cmake"

rm -rf ../cmake_build

cmake -B ../cmake_build -S ../. \
  -DUSE_MPI=ON \
  -DUSE_SUNDIALS=ON \
  -DUSE_MIZUROUTE=ON \
  -DSPECIFY_LAPACK_LINKS=OFF \
  -DCMAKE_BUILD_TYPE=Release

cmake --build ../cmake_build --target all -j 8
