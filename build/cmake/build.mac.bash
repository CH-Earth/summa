#!/bin/bash

# build SUMMA on a Mac using Bash, from cmake directory run this as ./build.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary




# Mac Example using homebrew:
export FC=/opt/homebrew/bin/gfortran                             # Fortran compiler family
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=../../../sundials/build/

cmake -B ../cmake_build -S ../. \
    -DUSE_MPI=ON \
    -DUSE_SUNDIALS=ON \
    -DUSE_MIZUROUTE=ON \
    -DSPECIFY_LAPACK_LINKS=ON \
    -DCMAKE_BUILD_TYPE=Release

cmake --build ../cmake_build --target all -j
