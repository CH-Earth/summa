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

If you want to use Sundials IDA or BE Kinsol, set -DUSE_SUNDIALS=ON in the build script.  Then, before summa can be built, Sundials needs to be installed. 

### Installing SUNDIALS
Download the file `sundials-X.Y.Z.tar.gz` (where X.Y.Z is the latest SUNDIALS version) at https://github.com/LLNL/sundials/releases/latest. Move `sundials-X.Y.Z.tar.gz` into `top_dir`. In other words, `$ ls top_dir` should include `summa sundials-X.Y.Z.tar.gz. Download this to the folder your preferred ${SUN_DIR}.

Extract the corresponding compressed file and rename
    $ cd sun_dir
    $ tar -xzf sundials-X.Y.Z.tar.gz && mv sundials-X.Y.Z sundials-software
    $ rm sundials-X.Y.Z.tar.gz

Create new empty directories to prep for SUNDIALS installation, within ${SUN_DIR}:
    $ mkdir sundials
    $ cd sundials
    $ mkdir builddir instdir

Copy CMake build script from SUMMA files to properly configure SUNDIALS from your chosen ${SUMMA_DIR}
    $ cd builddir
    $ cp ${SUMMA_DIR}/summa/build/cmake_external/build_cmakeSundials.bash .

Build SUNDIALS configured for SUMMA, within `builddir`: 
    $ ./build_cmakeSundials.bash
    $ make
    $ make install

We suggest you periodically update to the latest version. It is also possible to install using Git:  
    $ git clone https://github.com/LLNL/sundials.git sundials-software
    $ cd sundials-software
    $ git fetch --all --tags --prune
    $ git checkout tags/vX.Y.Z     
    
Note if you need to recompile after a system upgrade, delete the contents of sundials/instdir and sundials/buildir EXCEPT sundials/buildir/build_cmakeSundials.bash before building and installing.

### Building and installing SUMMA
First, you will need to tell CMake where Sundials is if you installed it and plan to use it:
    $ export CMAKE_PREFIX_PATH=/sun_dir/sundials/instdir

Then, a CMake build system must be generated.  E.g., from the top `summa/build/cmake` directory, using SUNDIALS:
    $ cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON

After there is build system directory, the shared library can be built using the `summabmi` CMake target. For example, the SummaSundials shared library file (i.e., the build config's `summabmi` target) can be built using:
    $ cmake --build ../cmake_build --target all -j

This will build a `cmake_build/libsumma.<version>.<ext>` file, where the version is configured within the CMake config, and the extension depends on the local machine's operating system.    

There is an example of a bash script to build the summa libraries at /build/cmake/build[_actors].[system_type].bash. SUNDIALS is turned on here. These need to be run in the cmake directory.

