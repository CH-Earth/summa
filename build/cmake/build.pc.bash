#!/bin/bash

# Build SUMMA on a PC using Bash, from cmake directory run this as ./build.pc.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# PC Example using Ubuntu: LAPACK Builds
#export FLAGS_OPT="-flto=1;-fuse-linker-plugin"                # optional compiler flags -- LAPACK builds

# PC Example using Ubuntu: Intel oneMKL builds (see https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)
#oneAPI_dir=/opt/intel/oneapi                                                 # Intel oneAPI main directory
#source $oneAPI_dir/setvars.sh                                                # initialize environment variables for Intel oneAPI
#export LIBRARY_LINKS="-m64;-Wl,--start-group ${MKLROOT}/lib/libmkl_gf_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group;-lpthread;-lm;-ldl" # static sequential library (i.e., no multithreading)
#export FLAGS_OPT="-m64;-I"${MKLROOT}/include";-flto=1;-fuse-linker-plugin"   # optional compiler flags -- Intel oneMKL builds

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DUSE_OPENWQ=OFF -DSPECIFY_LAPACK_LINKS=ON #-DCMAKE_BUILD_TYPE=Debug
cmake --build ../cmake_build --target all -j 
