#!/bin/bash
  
# build nextgen on Copernicus, from ngen directory put this one directory up and run this as ../build_ngen.cluster.bash
# for Summa
module load StdEnv/2023
module load gcc/12.3
module load openblas/0.3.24
module load openmpi/4.1.5
module load netcdf-fortran/4.6.1

export FLAGS_OPT="-flto=1;-fuse-linker-plugin"
export SUNDIALS_PATH=../../../sundials/instdir/

# for NextGen
module load boost
module load udunits/2.2.28

cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DUSE_NEXTGEN=ON
cmake --build extern/summa/cmake_build --target all -j

cmake -B cmake_build -S . -DNGEN_WITH_MPI:BOOL=OFF -DNGEN_WITH_PYTHON:BOOL=OFF -DNGEN_WITH_BMI_C:BOOL=ON -DNGEN_WITH_BMI_FORTRAN:BOOL=ON -DNGEN_WITH_NETCDF:BOOL=OFF
# can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb
# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, also turn MPI on
make -C cmake_build
