# SUMMA Installation

We have successfully installed SUMMA on a number of Unix-like (\*nix) operating systems, including Linux and Darwin (Mac OS X). Since we do a lot of our development on OS X, we have a [separate page](SUMMA_on_OS_X.md) on how to install the necessary tools and libraries on that platform. If you do not want to deal with installing programs and libraries and just want to run SUMMA, then we also have a SUMMA release that uses [Docker](https://www.docker.com). Details can be found on our [SUMMA using Docker](SUMMA_docker.md) page. If you plan to use Docker, then you can skip the rest of this page.

## Dependencies
To compile SUMMA, you will need (longer descriptions at the [bottom](#extended-description-of-dependencies) of this page):

* A Fortran compiler. The open source GNU Fortran compiler (`gfortran`) is a good choice.
* The [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) libraries, compiled/installed with Fortran support (i.e., `libnetcdff.*` is available on the system)
* The [LAPACK](http://www.netlib.org/lapack/) (Linear Algebra PACKage) library
* A copy of the SUMMA source code

Optional but highly recommended:

* The [SUNDIALS](https://sundials.readthedocs.io/en/latest/index.html) library, to benefit from recent advances to SUMMA's numerical implementation. Without the SUNDIALS library, the modeling options `num_method = kinsol` and `num_method = ida` won't work. SUNDIALS requires specific install instructions because by default its fortran module won't be build. A typical install will look as follows:

```bash
git clone https://github.com/LLNL/sundials.git sundials-software
cd sundials-software
git fetch --all --tags --prune
git checkout tags/vX.Y.Z (use most recent tag)
cd sundials
mkdir builddir && mkdir instdir && cd builddir
cmake ../../sundials-software/ -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples
make
make install
```

## Compilation
The preferred method to compile SUMMA is with `CMake`. This enables parallelization for faster builds and is generally easier than the [legacy approach using a Makefile](#legacy-makefile-instructions). CMake instructions can be found below.

### CMake

The SUMMA repository contains a number of CMake scripts in the `summa/build/cmake` directory. These scripts use CMake to compile SUMMA using the general `CMakeLists.txt` file found in the `summa/build` directory. The scripts compile SUMMA with a number of different options and for different OSX scenarios. Briefly:

- `FindNetCDF.cmake`: used to find the location (i.e., system path) for the Fortran NetCDF library. This path is needed in the actual compile scripts. 
- `FindOpenBLAS.cmake`: used to find the location (i.e., system path) for the OpenBLAS library. This path is needed in the actual compile scripts.
- `build.cluster.bash`:  compile SUMMA with SUNDIALS support on Digital Research Alliance Canada (DRAC) or similar infrastructure (e.g. Graham, Fir).
- `build.mac.bash`:  compile SUMMA with SUNDIALS support on macOS. Assumes MacPorts as local library manager (see example below).
- `build.pc.bash`:  compile SUMMA with SUNDIALS support on Windows. Experimental.
- `build_ngen.cluster.bash`: compile SUMMA with SUNDIALS support and NextGen integration on DRAC or similar. See specific instructions inside script.
- `build_ngen.mac.bash`: compile SUMMA with SUNDIALS support and NextGen integration on macOS. See specific instructions inside script.
- `summabmi.pc.in`: support file for NextGen integration.

### Submodules

Some build options need a submodule. A plain build needs none of them, so fetch only what you
turn on:

| Submodule | Needed by |
| --- | --- |
| `external/mizuRoute` | `USE_MIZUROUTE=ON` |
| `external/toml-f` | `USE_TOML=ON`, which `USE_MIZUROUTE` and `USE_MPI` turn on for you |

```bash
git submodule update --init external/mizuRoute external/toml-f
```

### Build options

Options default to `OFF`, except `USE_NGEN_FORCING` and `USE_NGEN_OUTPUT`, which default to `ON` and only take effect when using `USE_NEXTGEN=ON`. Set each option with `-DOPTION=ON` or `-DOPTION=OFF`.

| Option | Effect |
| --- | --- |
| `USE_SUNDIALS` | Build with the IDA and KINSOL solvers. Required for `num_method` of `ida` or `kinsol`. Needs `SUNDIALS_DIR`. |
| `USE_MPI` | Additionally build an MPI executable that distributes GRUs across ranks, and the calibration executable. |
| `USE_MIZUROUTE` | Build with mizuRoute river-network routing. |
| `USE_TOML` | Build the TOML configuration reader, needed for `-c`. Turned on automatically by `USE_MIZUROUTE` and `USE_MPI`. |
| `USE_NEXTGEN` | Build the BMI library for the NextGen framework instead of the executables. |
| `USE_NGEN_FORCING` | With `USE_NEXTGEN=ON`, use framework-supplied forcing. Set to `OFF` to use SUMMA's native NetCDF forcing reader. |
| `USE_NGEN_OUTPUT` | With `USE_NEXTGEN=ON`, use framework-managed output and skip SUMMA's native output writer. Set to `OFF` to retain native output. |
| `USE_OPENWQ` | Build with the OpenWQ water-quality coupler. |
| `SPECIFY_LAPACK_LINKS` | Take LAPACK flags from the `LIBRARY_LINKS` environment variable rather than detecting them. |

The executable name records the options, so differently configured builds share `bin/` without
overwriting one another: `summa.exe`, `summa_sundials.exe`, `summa_sundials_mizuroute.exe`, and so
on, plus `summa_sundials_mpi.exe` and `summa_sundials_opt.exe` from `USE_MPI`. `CMAKE_BUILD_TYPE`
is *not* part of the name, so a Debug and a Release build of the same options still collide; give
them separate build directories.

Each build script also takes a `clean` argument, which runs the CMake clean target and removes the
build directory:

```bash
./build.mac.bash clean
```

Building SUMMA-Actors is done from the [Summa-Actors](https://github.com/uofs-simlab/Summa-Actors)
repository, which carries its own build scripts and includes SUMMA as a submodule.

Most users will be able to compile SUMMA using one of the scripts above, after ensuring the paths in the scripts are set appropriately. As an example, imagine you're compiling on macOS but use Homebrew to manage your libraries, and that you installed SUNDIALS somewhere that's not the top-level summa folder. In this case, the default (MacPorts) path to the `gfortran` compiler, as well as the path to the SUNDIALS install directory, in the `build.mac.bash` script are not correct. You would need to update the script as follows:

```bash
#!/bin/bash
export FC=/opt/homebrew/bin/gfortran                             # Fortran compiler family
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=/your/path/to/sundials/instdir/

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DSPECIFY_LAPACK_LINKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build ../cmake_build --target all -j
```

You can test if SUMMA was compiled successfully by navigating to the new `bin` directory and running the newly created executable without any command line arguments. If the compilation is successful, you will see the help output as shown below:

```bash
> cd ~/path/to/summa/bin
> ./summa_sundials.exe

Usage: summa_sundials.exe [-m control_file] [-c config_file] [--manifest manifest_file] [-n newFileFreq] [-s fileSuffix] [-g startGRU countGRU] [-h iHRU] [-r freqRestart] [-p freqProgress] [--param name value]

Running options:
 -m --control       Define path/name of legacy SUMMA control file
 -c --config        Define path/name of TOML configuration file
 --manifest         Define path/name of multi-case run manifest
                     - At least one of --control, --config, or --manifest is required
 -n --newFile       Define frequency [noNewFiles,newFileEveryOct1] of new output files
 -s --suffix        Add fileSuffix to the output files
 -g --gru           Run a subset of countGRU GRUs starting from index startGRU
 -h --hru           Run a single HRU with index of iHRU
 -r --restart       Define frequency [y,m,d,e,never] to write restart files
 -p --progress      Define frequency [m,d,h,never] to print progress
 -v --version       Display version information of the current build
```

Continue reading [SUMMA configuration](../configuration/SUMMA_configuration.md) to learn more about how to configure SUMMA for your application. We strongly recommend that you get the [test applications](SUMMA_test_cases.md) to help you get started.

## Extra

### Extended description of dependencies
This is an extended description of the shorter list described above. To compile SUMMA, you will need:

 * a Fortran compiler. We have successfully used the intel Fortran compiler (`ifort`, version 17.x) and the GNU Fortran compiler (`gfortran`, version 6 or higher), the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile SUMMA with other Fortran compilers as well.

    If you do not have a Fortran compiler, you can install `gfortran` for free. The easiest way is to use a package manager. Note that `gfortran` is installed as part of the `gcc` compiler suite.

 * the NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. They are widely used in the hydrometeorological community and eventually almost all SUMMA I/O will use NetCDF. Most \*nix package managers include a NetCDF port. Note that you need to ensure that:

    * You have NetCDF version 4.x;
    * The NetCDF libraries are compiled with the same compiler as you plan to use for compiling SUMMA; and
    * You have the NetCDF Fortran library installed (`libnetcdff.*`) and not just the C-version.

 * the [LAPACK](http://www.netlib.org/lapack/) (Linear Algebra PACKage) library provides a series of routines for linear algebra operations, including matrix solvers. How to install the library depends on your \*nix variant and is not covered here. For example, on OS X you will get all the necessary LAPACK routines by installing the ATLAS software (again, this is easiest using a package manager; note that ATLAS can take many hours to build the first time when you install it).

 * a copy of the SUMMA source code from [this repo](https://github.com/NCAR/summa). You have a number of options:

    * If you just want to use the latest stable release of SUMMA, then simply look for the [latest release](https://github.com/NCAR/summa/releases);
    * If you want the latest and greatest (and potentially erroneous), download a copy of the [development branch](https://github.com/ncar/summa/tree/develop) (or clone it);
    * If you may want to do SUMMA development, then fork the repo on github and start editing your own copy.
    
    
### Legacy Makefile instructions
These instructions have been included here for legacy users. If possible, use the CMake instructions for an easier install experience. 

Once you have all the listed dependencies, you can compile SUMMA using the following steps for using the `Makefile`:

 1. Navigate to your local copy of the SUMMA directory and go to the `build` subdirectory;

 1. Specify a number of environment variables that are used by the build process. You will need to set the following:

    * `F_MASTER` : This is the parent directory of the `build` directory.
    
        > Example: Given the following directory structure: 
        >            ```
        >
        >            summa/  
        >            ├── build  
        >            ├── ci  
        >            ├── COPYING  
        >            ├── Dockerfile  
        >            ├── docs  
        >            ├── header.license  
        >            ├── LICENSE.txt  
        >            ├── mkdocs.yml  
        >            ├── readme.md -> docs/index.md  
        >            └── utils  
        >            ```
        >
        >  `export F_MASTER=/summa` 
        
    * `FC` : Define the compiler family. This is only used to determine the compiler flags. If your compiler is not included, you will need to add the relevant section to the `Makefile`.
    
        > Example: `export FC=gfortran`
            
    * `FC_EXE` : This is the actual compiler executable that is invoked.
    
        > Example: `export FC_EXE=gfortran`
        
    * `INCLUDES`: This is the path to the NetCDF and LAPACK include files. This is typically `/usr/include` or `/usr/local/include`.   
        
        > Example: `export INCLUDES='-I/usr/include -I/usr/local/inclde -I<other paths>`
            
    * `LIBRARIES`: This is the path to the NetCDF and LAPACK libraries, typically `/usr/lib`. The following command will help you determine the correct paths: `find / -type f \( -name "libblas*.so*" -o -name "libnetcdf*.so*" \)`. 
        
        > Example: `export LIBRARIES='-L/usr/lib -lnetcdff -L/usr/lib/x86_64-linux-gnu -llapack -lblas'`

    If you are using the `bash` shell, then you would set these environment variables with `export FC=gfortran` for example. You may need to modify the `Makefile` if you are using a different Fortran compiler or your setup is different. If someone wants to contribute an actual `configure` script that would be great.

    * If you are compiling SUMMA using packages installed through `Homebrew` with `gfortran`, then use the following entries in Part 0 of Makefile:
    Date updated: May-2024
        ```
        FC = gfortran
        FC_EXE = gfortran
        INCLUDES = -I/opt/homebrew/Cellar/netcdf-fortran/x.x.x/include -I/opt/homebrew/Cellar/lapack/x.x.x/include
        LIBRARIES = -L/opt/homebrew/Cellar/netcdf-fortran/x.x.x/lib -lnetcdff -L/opt/homebrew/Cellar/lapack/x.x.x/lib -lblas -llapack
        ```
    > Note: change `x.x.x` with the exact version number in both the `INCLUDE` and `LIBRARIES` variables. 
    This can be done by using `ls /opt/homebrew/Cellar/netcdf-fortran/` and then using `tab` button to find the current version installed on your machine. 
    Do the same for `lapack`. Currently, the most up-to-date `netcdf-fortran` version is `4.6.1` and lapack is `3.12.0`.

    * If you are compiling SUMMA on the [Graham cluster](https://docs.alliancecan.ca/wiki/Graham/en) of the [Digital Research Alliance of Canada](https://alliancecan.ca/en) using `ifort` then use the following entries in Part 0 of Makefile:
    Date updated:July-2020
        ```
        FC = ifort
        FC_EXE = /cvmfs/restricted.computecanada.ca/easybuild/software/2017/Core/ifort/2019.3.199/compilers_and_libraries_2019.3.199/linux/bin/intel64/ifort
        INCLUDES = -I/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx/Compiler/intel2016.4/netcdf-fortran/4.4.4/include
        LIBRARIES = -L/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx/Compiler/intel2016.4/netcdf-fortran/4.4.4/lib -lnetcdff -mkl
        ```
    > Note: Before compiling, load the most recent intel module. Check the available versions with: 
`module spider intel`. Then load the latest version: `module load intel/2019.3`. 
Lapack and blas libraries are loaded with library argument `-mkl`. 
------------------------------



 1. Check that all variables in the Makefile are set correctly by typing `make check`. Inspect the variables and make sure that they make sense. If not, modify the Makefile further.

 1. Type `make` (if you are in the `build` directory). If all goes well, this will build SUMMA and move the executable `summa.exe` to the `bin` directory. You may get some warnings (depending on your compiler settings), but you should not get any errors;

 1. Pay attention to the `make` output. You may need to set some environment variables (`LD_LIBRARY_PATH` in particular) to support dynamic linking;

 1. If the code compiles successfully, then the last line of output from the make process will tell you where the SUMMA executable is installed (it goes into `summa/bin`). Run `summa.exe` in that directory (you may need to provide the full path). If all goes well, you should get a usage message that looks something like (depending on the SUMMA version):

```
    Usage: summa.exe -m master_file [-s fileSuffix] [-g startGRU countGRU] [-h iHRU] [-r freqRestart] [-p freqProgress] [-c]
     summa.exe          summa executable

    Running options:
     -m --master        Define path/name of master file (required)
     -s --suffix        Add fileSuffix to the output files
     -g --gru           Run a subset of countGRU GRUs starting from index startGRU
     -h --hru           Run a single HRU with index of iHRU
     -r --restart       Define frequency [y,m,d,never] to write restart files
     -p --progress      Define frequency [m,d,h,never] to print progress
     -v --version       Display version information of the current built
```

If you get this far then SUMMA is installed correctly and functional.
