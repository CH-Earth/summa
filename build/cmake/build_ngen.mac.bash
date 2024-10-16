#!/bin/bash
  
# build nextgen on Mac, from ngen directory put this one directory up and run this as ../build_ngen.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use
export C_INCLUDE_PATH=/opt/local/include
export CPLUS_INCLUDE_PATH=/opt/local/include
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=../../../sundials/instdir/

cmake -B extern/summa/cmake_build -S extern/summa -DUSE_NEXTGEN=ON -DSPECIFY_LAPACK_LINKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build extern/summa/cmake_build --target all -j

cmake -B cmake_build -S . -DBoost_INCLUDE_DIR=/opt/local/libexec/boost/1.81/include -DNGEN_WITH_BMI_FORTRAN=ON -DPython_NumPy_INCLUDE_DIR=/opt/local/bin/python -DNGEN_WITH_MPI:BOOL=OFF -DNGEN_WITH_PYTHON:BOOL=OFF -DNGEN_WITH_BMI_C:BOOL=ON -DNGEN_WITH_BMI_FORTRAN:BOOL=ON -DNGEN_WITH_NETCDF:BOOL=OFF
# can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb
# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, also turn MPI on
make -C cmake_build
