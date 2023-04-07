#! /bin/bash

#### parent directory of the 'build' directory ####
export ROOT_DIR=/Users/amedin/Research/SummaBMI/ngen/extern/summa

export INCLUDES="-I/opt/local/include -I/opt/local/lib"
export LDFLAGS="-L/opt/local/lib"
export LIBRARIES="${LDFLAGS} -llapack -lgfortran -lnetcdff -lnetcdf"

make -f ${ROOT_DIR}/summa/build/makefiles/makefile_mac
export LD_LIBRARY_PATH=${ROOT_DIR}/summa/bin
