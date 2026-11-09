#!/bin/bash

# from {$maindir}/sundials/builddir, run
# cp ../../summa/build/build_cmakeSundials.bash build_cmake.bash
# run script from the builddir directory with ./build_cmake.bash
# Note, using -DSUNDIALS_ENABLE_C_EXAMPLES=OFF -DSUNDIALS_ENABLE_FORTRAN_EXAMPLES=OFF for quickness

cmake ../../sundials-software/ -DSUNDIALS_ENABLE_C_EXAMPLES=OFF -DSUNDIALS_ENABLE_FORTRAN_EXAMPLES=OFF -DSUNDIALS_ENABLE_FORTRAN=ON -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples
