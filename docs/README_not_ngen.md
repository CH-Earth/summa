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

### Solver options in SUMMA build scripts
First, cd into the cmake directory in summa, without Actors:

    cd summa/build/cmake

or with Actors:

    cd summa/build/summa/build/cmake

If you want to use Sundials IDA or BE Kinsol, set -DCMAKE_BUILD_TYPE=Sundials* in the build script instead of BE*.  Then, before summa can be built, Sundials needs to be installed.

### Installing SUNDIALS
0. For reference, let the main summa directory be contained in the folder `top_dir`
    - e.g., summa exectuables would be located within `top_dir/summa/bin`
    
2. Download the file `sundials-X.Y.Z.tar.gz` (where X.Y.Z is the latest SUNDIALS version) at https://github.com/LLNL/sundials/releases/latest 
 
3. Move `sundials-X.Y.Z.tar.gz` into `top_dir`
    - `$ ls top_dir` should include `summa sundials-X.Y.Z.tar.gz` 

4. Extract the corresponding compressed file and rename
    1. `$ cd top_dir`
    2. `$ tar -xzf sundials-X.Y.Z.tar.gz && mv sundials-X.Y.Z sundials-software`
    3. `sundials-X.Y.Z.tar.gz` can now be removed to save space if desired

5. Create new empty directories to prep for SUNDIALS installation
    1. within `top_dir`: `$ mkdir sundials`
    2. `$ cd sundials`
    3. `$ mkdir builddir instdir`

6. Copy CMake build script from SUMMA files to properly configure SUNDIALS
    1. `$ cd builddir`
    2. `$ cp ../../summa/build/cmake_external/build_cmakeSundials.bash .`

7. Build SUNDIALS configured for SUMMA
    1. within `builddir`: `$ ./build_cmakeSundials.bash`
    2. `$ make`
    3. `$ make install`


We suggest you periodically update to the latest version. It is also possible to install using Git: 
 
    $ git clone https://github.com/LLNL/sundials.git sundials-software
    $ cd sundials-software
    $ git fetch --all --tags --prune
    $ git checkout tags/vX.Y.Z     
    
Note if you need to recompile after a system upgrade, delete the contents of sundials/instdir and sundials/buildir EXCEPT sundials/buildir/build_cmake before building and installing.

Note that when there is an existing directory, it may sometimes be necessary to clear it and regenerate, especially if any changes were made to the CMakeLists.txt file.

### Building SUMMA

After there is build system directory, the shared library can be built using the `summabmi` CMake target. For example, the SummaSundials shared library file (i.e., the build config's `summabmi` target) can be built using:

    $ cmake --build ../cmake_build --target all

This will build a `cmake_build/libsummabmi.<version>.<ext>` file, where the version is configured within the CMake config, and the extension depends on the local machine's operating system.    

There is an example of a bash script to build the summa libraries at /build/cmake/build[_actors].[system_type].bash. Sundials is turned on here. These need to be run in the cmake directory.

