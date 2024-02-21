#!/bin/bash
  
# build on Mac, from cmake directory run this as ./build_actors.mac.bash

# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
export LINK_DIRS=/opt/local/lib                               # Link directories for cmake
export INCLUDES_DIRS='/opt/local/include;/opt/local/lib'      # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
export LIBRARY_LINKS='-llapack;-lgfortran;-lnetcdff;-lnetcdf' # list of library links (cmake uses semicolons as separators)
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials_Actors
cmake --build ../cmake_build --target all
