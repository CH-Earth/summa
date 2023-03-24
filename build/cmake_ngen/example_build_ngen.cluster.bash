#!/bin/bash
  
# build nextgen on Copernicus, from ngen directory put this one directory up and run this as ../example_build_ngen.bash

#module load gnu
module load boost
module load udunits/2.2.28
#module load python/3.7.5
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2


cmake -B /extern/cfe/cmake_build -S /extern/cfe
cmake --build /extern/cfe/cmake_build --target all

cmake -B /extern/topmodel/cmake_build -S /extern/topmodel
cmake --build /extern/topmodel/cmake_build --target all

cmake -B /extern/sloth/cmake_build -S /extern/sloth
cmake --build /extern/sloth/cmake_build --target all

cmake -B /extern/noah-owp-modular/cmake_build -S /extern/noah-owp-modular
cmake --build /extern/noah-owp-modular/cmake_build --target all

cmake -B /extern/iso_c_fortran_bmi/cmake_build -S /extern/iso_c_fortran_bmi
cmake --build /extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B /extern/snow17/cmake_build -S /extern/snow17
cmake --build /extern/snow17/cmake_build --target all

cmake -B /extern/summa/cmake_build -S /extern/summa
cmake --build /extern/summa/cmake_build --target all

cmake -B /extern/evapotranspiration/cmake_build -S /extern/evapotranspiration
cmake --build /extern/evapotranspiration/cmake_build --target all

cmake -B cmake_build -S . -DMPI_ACTIVE:BOOL=ON -DNGEN_ACTIVATE_PYTHON:BOOL=OFF -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON  \
#cmake -B cmake_build -S .  -DNGEN_ACTIVATE_PYTHON:BOOL=OFF -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON  \
#  -DUDUNITS2_LIBRARY=/glade/u/apps/ch/opt/udunits/2.2.28/gnu/10.1.0/lib/libudunits2.so \
#  -DUDUNITS2_INCLUDE=/glade/u/apps/ch/opt/udunits/2.2.28/gnu/10.1.0/include # \
#  -DCMAKE_BUILD_TYPE=Debug

  # can also add -DCMAKE_BUILD_TYPE=Debug to be able to run in gdb

#make -j 8 -C cmake_build    # build w/ 8 parallel jobs
make -C cmake_build

