#!/bin/bash
  
# build on Copernicus or Graham, from cmake directory run this as ./build.cluster.bash
# for Summa
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials_Cluster
cmake --build ../cmake_build --target all

