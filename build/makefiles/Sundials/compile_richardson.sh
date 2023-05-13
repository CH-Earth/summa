#! /bin/bash

#### parent directory of the 'build' directory ####
export ROOT_DIR=/u1/gwu479/SummaSundials

export INCLUDES="-I/usr/include -I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"
export LIBRARIES="${LDFLAGS} -lnetcdff -lopenblas -lnetcdf"

make -f ${ROOT_DIR}/summa/build/makefiles/makefile_cluster
export LD_LIBRARY_PATH=${ROOT_DIR}/summa/bin
