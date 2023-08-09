#!/bin/bash
  
# build on Mac, from cmake directory run this as ./build_actors.mac.bash
export FC=gfortran                                            # Fortran compiler family
export LINK_DIRS='/usr/lib;/usr/lib/x86_64-linux-gnu'                              # Link directories for cmake
export INCLUDES_DIRS='/usr/include;/usr/local/include'      # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
export LIBRARY_LINKS='-llapack;-lgfortran;-lnetcdff;-lnetcdf' # list of library links (cmake uses semicolons as separators)
export DIR_SUNDIALS='/usr/local/sundials-6.3.0'
cmake -B ../cmake_build_regular -S . -DCMAKE_BUILD_TYPE=Sundials_Actors
cmake --build ../cmake_build_regular --target all

