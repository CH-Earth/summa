## Compile SUMMA -- assumes that this script is in the build directory

## manually set environment variables
export FC=gfortran                                                                  # Fortran compiler family
export FC_EXE=gfortran # name of Fortran executable
export INCLUDES='-I/usr/include -I/usr/local/include'                               # path to NetCDF and LAPACK include files
export LIBRARIES='-L/usr/lib -lnetcdff -L/usr/lib/x86_64-linux-gnu -llapack -lblas' # path to libraries

## automatically set environment variables
export F_MASTER=$(dirname $PWD)                                                     # parent of build directory

## link and compile
make
