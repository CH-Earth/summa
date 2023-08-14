#!/bin/bash

# Build SUMMA on a PC using Bash, from cmake directory run this as ./build.pc.bash
# Environment variables that must be set for cmake: FC, LINK_DIRS, INCLUDES_DIRS, and LIBRARY_LINKS
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# PC Example using Ubuntu: LAPACK Builds
#export FC=gfortran                                            # Fortran compiler family
#export LINK_DIRS=/usr/local/lib                               # Link directories for cmake
#export INCLUDES_DIRS="/usr/include;/usr/local/include"        # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
#export LIBRARY_LINKS="-llapack;-lnetcdff;-lnetcdf"            # list of library links -- LAPACK builds
#export FLAGS_OPT="-fuse-linker-plugin"                                           # optional compiler flags -- LAPACK builds

# PC Example using Ubuntu: Intel oneMKL builds (see https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)
#oneAPI_dir=/opt/intel/oneapi                                  # Intel oneAPI main directory
#source $oneAPI_dir/setvars.sh                                 # initialize environment variables for Intel oneAPI
#export FC=gfortran                                            # Fortran compiler family
#export LINK_DIRS=/usr/local/lib                               # Link directories for cmake
#export INCLUDES_DIRS="/usr/include;/usr/local/include"        # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
#export LIBRARY_LINKS="-lnetcdff;-lnetcdf;-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group;-lgomp;-lpthread;-lm;-ldl"         # list of library links -- Intel oneMKL builds
#export FLAGS_OPT="--m64;-I"${MKLROOT}/include;-fuse-linker-plugin"                # optional compiler flags -- Intel oneMKL builds

# CMake Commands (build type controlled using DCMAKE_BUILD_TYPE)
cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=BE
cmake --build ../cmake_build --target all

