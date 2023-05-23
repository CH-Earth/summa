#!/bin/bash
  
# build nextgen on Copernicus, from ngen directory put this one directory up and run this as ../build_ngen.cluster.bash
# for Summa
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2

# for NexGen
module load boost
module load udunits/2.2.28

cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DCMAKE_BUILD_TYPE=BE_NexGen_Cluster
cmake --build extern/summa/cmake_build --target all

cmake -B cmake_build -S . -DMPI_ACTIVE:BOOL=OFF -DNGEN_ACTIVATE_PYTHON:BOOL=OFF -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON
# can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb
# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, also turn MPI on
make -C cmake_build

