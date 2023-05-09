#!/bin/bash
  
# build on Mac, from cmake directory run this as ./build.mac.bash

cmake -B ../cmake_build -S summa -DCMAKE_BUILD_TYPE=IDA
cmake -j ../cmake_build --target all

