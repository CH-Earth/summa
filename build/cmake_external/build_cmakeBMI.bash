#!/bin/bash

# from {$maindir}/bmi/builddir, run
# cp ../../summa/build/cmake_external/build_cmakeBMI.bash build_cmake
# run script from the builddir directory with ./build_cmake
# run `make`, then `make install`

cmake ../../bmi-fortran/ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../bmi/instdir
