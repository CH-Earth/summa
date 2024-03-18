#!/bin/bash
  
# build nextgen on Mac, from ngen directory put this one directory up and run this as ../build_ngen.mac.bash
# Mac Example using MacPorts:
export FC=/opt/local/bin/gfortran                             # Fortran compiler family
export LINK_DIRS=/opt/local/lib                               # Link directories for cmake
export INCLUDES_DIRS='/opt/local/include;/opt/local/lib'      # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
export LIBRARY_LINKS='-llapack;-lgfortran;-lnetcdff;-lnetcdf' # list of library links (cmake uses semicolons as separators)
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use

cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DCMAKE_BUILD_TYPE=BE_NexGen
cmake --build extern/summa/cmake_build --target all

cmake -B cmake_build -S . -DBoost_INCLUDE_DIR=/opt/local/libexec/boost/1.81/include -DPython_NumPy_INCLUDE_DIR=/opt/local/bin/python -DNGEN_WITH_MPI:BOOL=OFF -DNGEN_WITH_PYTHON:BOOL=OFF -DNGEN_WITH_BMI_C:BOOL=ON -DNGEN_WITH_BMI_FORTRAN:BOOL=ON -DNGEN_WITH_NETCDF:BOOL=OFF
# can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb
# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, also turn MPI on
make -C cmake_build
