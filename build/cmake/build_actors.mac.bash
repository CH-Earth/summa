#!/bin/bash
  
# build on Mac, from cmake directory run this as ./build_actors.mac.bash

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials_Actors
cmake --build ../cmake_build --target all

