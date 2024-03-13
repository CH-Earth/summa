# SummaSundials 

## About

This is configured with the [CMakeLists.txt](CMakeLists.txt) and other files in this directory.

#### Extra Outer Directory for Actors

Actors needs one


### Getting the Latest Changes

There are two steps to getting upstream submodule changes fully 
  1. fetching and locally checking out the changes from the remote
  2. committing the new checkout revision for the submodule

To fetch and check out the latest revision (for the [currently used branch](#viewing-the-current-branch)):

    git pull

# Usage

## Building Libraries

First, cd into the cmake directory in summa, without Actors:

    cd summa/build/cmake

or with Actors:

    cd summa/build/summa/build/cmake

If you want to use Sundials IDA or BE Kinsol, set -DCMAKE_BUILD_TYPE=Sundials* in the build script instead of BE*.  Then, before summa can be built, Sundials needs to be installed. 
Download the latest release of IDA solver from SUNDIALS package in https://computing.llnl.gov/projects/sundials/sundials-software, using X.Y.Z as the latest version

    wget "https://github.com/LLNL/sundials/releases/download/vX.Y.Z/sundials-X.Y.Z.tar.gz"

Extract the corresponding compressed file, rename

    tar -xzf sundials-X.Y.Z.tar.gz && mv sundials-X.Y.Z sundials-software

We suggest you periodically update to the latest version-- you can also install through github 
    git clone https://github.com/LLNL/sundials.git sundials-software
    cd sundials-software
    git fetch --all --tags --prune
    git checkout tags/vX.Y.Z
    
An example build_cmake file is at summa/build/cmake_external/build_cmakeSundials.bash which you can copy to builddir as build_cmake. Then, enter the buildir and run

    cd sundials/buildir
    ./build_cmake
    make
    make install
    
Note if you need to recompile after a system upgrade, delete the contents of sundials/instdir and sundials/buildir EXCEPT sundials/buildir/build_cmake before building and installing.


Note that when there is an existing directory, it may sometimes be necessary to clear it and regenerate, especially if any changes were made to the CMakeLists.txt file.

After there is build system directory, the shared library can be built using the `summabmi` CMake target. For example, the SummaSundials shared library file (i.e., the build config's `summabmi` target) can be built using:

    cmake --build ../cmake_build --target all

This will build a `cmake_build/libsummabmi.<version>.<ext>` file, where the version is configured within the CMake config, and the extension depends on the local machine's operating system.    

There is an example of a bash script to build the summa libraries at /build/cmake/build[_actors].[system_type].bash. Sundials is turned on here. These need to be run in the cmake directory.

