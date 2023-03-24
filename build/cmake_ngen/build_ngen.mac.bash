#!/bin/bash
  
# build nextgen on Mac, from ngen directory put this one directory up and run this as ../example_build_ngen.mac.bash

cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa
cmake --build extern/summa/cmake_build --target all

cmake -B /ngen/cmake_build -S . -DMPI_ACTIVE:BOOL=ON -DNGEN_ACTIVATE_PYTHON:BOOL=OFF -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON  \
#cmake -B /ngen/cmake_build -S .  -DNGEN_ACTIVATE_PYTHON:BOOL=OFF -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON  \
#  -DUDUNITS2_LIBRARY=/glade/u/apps/ch/opt/udunits/2.2.28/gnu/10.1.0/lib/libudunits2.so \
#  -DUDUNITS2_INCLUDE=/glade/u/apps/ch/opt/udunits/2.2.28/gnu/10.1.0/include # \
#  -DCMAKE_BUILD_TYPE=Debug

  # can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb

#make -j 8 -C cmake_build    # build w/ 8 parallel jobs
make -C cmake_build

