#!/bin/bash
  
# build on HPC, from cmake directory run this as ./build.cluster.bash

# load these modules in run environment as well as build environment
# Digital Resource Alliance of Canada settings
module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1
#
# Purdue Anvil settings
#module load gcc/14.2.0
#module load openmpi/4.1.6
#module load openblas/0.3.17
#module load netcdf-fortran/4.5.3

# Actors install of sundials
#export SUNDIALS_DIR="$CMAKE_PREFIX_PATH:$HOME/Summa-Actors/utils/dependencies/install/sundials/"
#
# Regular install of sundials
export SUNDIALS_DIR=$HOME/SummaSundials/sundials/instdir/

# May want to use this flag
#export FLAGS_OPT="-flto=1;-fuse-linker-plugin"

cmake -B ../cmake_build -S ../. \
    -DUSE_SUNDIALS=ON \
    -DUSE_MPI=OFF \
    -DUSE_NEXTGEN=OFF \
    -DUSE_OPENWQ=OFF \
    -DUSE_MIZUROUTE=OFF \
    -DSPECIFY_LAPACK_LINKS=OFF \
    -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake_build --target all -j
