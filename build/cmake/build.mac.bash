#!/bin/bash

# build SUMMA on a Mac using Bash, from cmake directory run this as ./build.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use

export SUNDIALS_DIR=../../../sundials/instdir/
export LINK_DIRS=/opt/local/lib                               # Link directories for cmake
export INCLUDES_DIRS='/opt/local/include;/opt/local/lib'      # directories for INCLUDES \
export LIBRARY_LINKS='-llapack'                               # list of library links 

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON
cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake_build --target all -j
