#!/bin/bash

# Script to compile the SUMMA model using cmake
# This will call the CMakeLists.txt file in one directory above the current
# directory, and create a directory cmake_build where all of the build files
# will be stored. This script will build Summa in parallel.

cmake -B cmake_build -S ../. 
cmake --build cmake_build --target all -j 