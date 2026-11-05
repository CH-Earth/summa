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
# One GRU thread per core the rank was given; BLAS stays on the calling thread
# ----------------------------------------------------------------------

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=256M
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
