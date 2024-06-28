#!/bin/bash

# build SUMMA on a Mac using Bash, from cmake directory run this as ./build.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON
cmake --build ../cmake_build --target all -j
