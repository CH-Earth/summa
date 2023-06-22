#! /bin/bash

#### load modules if using Compute Canada or Copernicus ####
module load gcc/9.3.0
module load netcdf-fortran
module load openblas
module load caf


#### parent directory of the 'build' directory ####
# export ROOT_DIR=/home/kklenk/Summa-Projects/Summa-Actors

# export INCLUDES="-I$EBROOTNETCDFMINFORTRAN/include"
# export LIBRARIES="-L$EBROOTNETCDFMINFORTRAN/lib64\
#     -L$EBROOTOPENBLAS/lib\
#     -lnetcdff -lopenblas"

# # INCLUDES FOR Actors Component
# export ACTORS_INCLUDES="-I$EBROOTCAF/include\
#     -I$EBROOTNETCDFMINFORTRAN/include"

# export ACTORS_LIBRARIES="-L$EBROOTCAF/lib\
#     -L$EBROOTCAF/lib64\
#     -L$EBROOTNETCDFMINFORTRAN/lib64\
#     -L$EBROOTOPENBLAS/lib\
#     -L$ROOT_DIR/bin\
#     -lcaf_core -lcaf_io -lsumma -lopenblas -lnetcdff"

# make -f ${ROOT_DIR}/build/makefiles/makefile_cluster cpp
# export LD_LIBRARY_PATH=${ROOT_DIR}/bin