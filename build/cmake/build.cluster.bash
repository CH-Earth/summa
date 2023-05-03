#!/bin/bash
  
# build on Copernicus or Graham, from summa directory put this one directory up and run this as ../build.cluster.bash

module load boost
module load udunits/2.2.28
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2


cmake -B cmake_build -S summa -DCMAKE_BUILD_TYPE=Cluster -DSUNDIALS_ACTIVE:BOOL=ON -DACTORS_ACTIVE:BOOL=OFF
cmake --build summa/cmake_build --target all

make -C cmake_build

