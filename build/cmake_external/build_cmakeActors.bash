#!/bin/bash

# from {$maindir}/actor-framework, run
# cp ../../summa/build/summa/build/cmake_external/build_cmakeActors.bash build_cmake
# run script from the actor-framework directory with ./build_cmake
# run `cd build`, `make`, then `make install`

export CXX="g++"
./configure --prefix=../install
