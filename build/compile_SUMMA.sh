## Compile SUMMA -- Not fully automatic -- some variables must be manually set (see below)
#Prerequisites: 1) Linux terminal (tested on Ubuntu 22.04 LTS but may work on other Linux platforms)
#               2) This script must be in the build directory (e.g. ~/summa/build)
#How to use:    1) Manually set FC, FC_EXE, INCLUDES, and LIBRARIES
#               2) Save file and enter `source compile_SUMMA.sh` in the terminal

## manually set environment variables
export FC=gfortran                                                                  # Fortran compiler family
export FC_EXE=gfortran                                                              # name of Fortran executable
export INCLUDES='-I/usr/include -I/usr/local/include'                               # path to NetCDF and LAPACK include files
export LIBRARIES='-L/usr/lib -lnetcdff -L/usr/lib/x86_64-linux-gnu -llapack -lblas' # path to libraries

# ------------------------------- Automatic Past This Point ------------------------------------- 

## automatically set environment variables
export F_MASTER=$(dirname $PWD)                                                     # parent of build directory

## link and compile
make
