#! /bin/bash

#### load modules if using Compute Canada or Copernicus ####
module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2

#### parent directory of the 'build' directory ####
export ROOT_DIR=/project/gwf/gwf_cmt/ngen/ngen_code/ngen/extern/summa

export INCLUDES="-I${EBROOTNETCDFMINFORTRAN}/include"
export LDFLAGS="-L${EBROOTNETCDFMINFORTRAN}/lib64 -L${EBROOTOPENBLAS}/lib"
export LIBRARIES="${LDFLAGS} -lnetcdff -lopenblas"

make -f ${ROOT_DIR}/summa/build/makefiles/makefile_cluster
export LD_LIBRARY_PATH=${ROOT_DIR}/summa/bin
