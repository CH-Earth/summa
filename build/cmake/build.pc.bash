#!/bin/bash

# build SUMMA on a PC using Bash, from cmake directory run this as ./build.pc.bash
# Environment variables that must be set for cmake: FC, LINK_DIRS, INCLUDES_DIRS, and LIBRARY_LINKS

# PC Example using Ubuntu
#export FC=gfortran                                            # Fortran compiler family
#export LINK_DIRS=/usr/local/lib                               # Link directories for cmake
#export INCLUDES_DIRS='/usr/include;/usr/local/include'        # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
#export LIBRARY_LINKS='-llapack;-lgfortran;-lnetcdff;-lnetcdf' # list of library links (cmake uses semicolons as separators)

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=BE
cmake --build ../cmake_build --target all

