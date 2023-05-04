#!/bin/bash
  
# build on Mac, from summa directory put this one directory up and run this as ../build.mac.bash

cmake -B cmake_build -S summa -DCMAKE_BUILD_TYPE=IDA 
cmake --build summa/cmake_build --target all

make -C cmake_build
