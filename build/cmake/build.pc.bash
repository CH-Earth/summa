#!/bin/bash

cmake -B ../cmake_build -S ../. 
cmake --build ../cmake_build --target all -j 
