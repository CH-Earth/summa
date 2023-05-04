#!/bin/bash
  
# build on Copernicus or Graham, from summa directory put this one directory up and run this as ../build.cluster.bash

module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2

# for NexGen
module load boost
module load udunits/2.2.28

# for Actors
module load caf


cmake -B ../cmake_build -S ../../ -DCMAKE_BUILD_TYPE=IDA_Cluster
cmake --build ../cmake_build --target all

make -C cmake_build

