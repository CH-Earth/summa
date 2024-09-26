#!/bin/bash
  
# build on Copernicus or Graham, from cmake directory run this as ./build.cluster.bash

module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1

export FLAGS_OPT="-flto=1;-fuse-linker-plugin"
export SUNDIALS_PATH=../../../sundials/instdir/

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake_build --target all -j
