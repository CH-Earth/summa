#!/bin/bash

#!/bin/bash

# Mac example using Homebrew
export FC=/opt/homebrew/bin/gfortran
export LIBRARY_LINKS='-llapack'
export SUNDIALS_DIR=../../../sundials/build/

cmake -B ../cmake_build -S ../. \
    -DUSE_MPI=ON \
    -DUSE_SUNDIALS=ON \
    -DUSE_MIZUROUTE=ON \
    -DSPECIFY_LAPACK_LINKS=ON \
    -DCMAKE_BUILD_TYPE=Release

cmake --build ../cmake_build --target all -j
