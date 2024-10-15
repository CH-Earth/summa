#!/bin/bash
  
# build on Copernicus or Graham, from cmake directory run this as ./build_actors.cluster.bash
# for Summa
module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1

# for Actors
module load caf

export FLAGS_OPT="-flto=1;-fuse-linker-plugin"
export SUNDIALS_PATH=/globalhome/kck540/HPC/Libraries/sundials/v7.0/instdir

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DUSE_ACTORS=ON -DSPECIFY_LAPACK_LINKS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake_build --target all -j
