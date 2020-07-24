# SUMMA Installation

We have successfully installed SUMMA on a number of Unix-like (\*nix) operating systems, including Linux and Darwin (Mac OS X). Since we do a lot of our development on OS X, we have a [separate page](SUMMA_on_OS_X.md) on how to install the necessary tools and libraries on that platform. If you do not want to deal with installing programs and libraries and just want to run SUMMA, then we also have a SUMMA release that uses [Docker](https://www.docker.com). Details can be found on our [SUMMA using Docker](SUMMA_docker.md) page. If you plan to use Docker, then you can skip the rest of this page.

To compile SUMMA, you will need:

 * a Fortran compiler. We have successfully used the intel Fortran compiler (`ifort`, version 17.x) and the GNU Fortran compiler (`gfortran`, version 6 or higher), the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile SUMMA with other Fortran compilers as well.

    If you do not have a Fortran compiler, you can install `gfortran` for free. The easiest way is to use a package manager. Note that `gfortran` is installed as part of the `gcc` compiler suite.

 * the NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. They are widely used in the hydrometeorological community and eventually almost all SUMMA I/O will use NetCDF. Most \*nix package managers include a NetCDF port. Note that you need to ensure that:

    * You have NetCDF version 4.x;
    * The NetCDF libraries are compiled with the same compiler as you plan to use for compiling SUMMA; and
    * You have the NetCDF Fortran library installed (`libnetcdff.*`) and not just the C-version.

 * the [LAPACK](http://www.netlib.org/lapack/) (Linear Algebra PACKage) library provides a series of routines for linear algebra operations, including matrix solvers. How to install the library depends on your \*nix variant and is not covered here. For example, on OS X you will get all the necessary LAPACK routines by installing the ATLAS software (again, this is easiest using a package manager; note that ATLAS can take many hours to build the first time when you install it).

 * the [ATLAS](http://math-atlas.sourceforge.net/) (Automatically Tuned Linear Algebra Software) library. Note that this is required on OS X using the gfortran compiler to be able to use LAPACK. It's not clear that this is used on other \*nix machines.

 * a copy of the SUMMA source code from [this repo](https://github.com/NCAR/summa). You have a number of options:

    * If you just want to use the latest stable release of SUMMA, then simply look for the [latest release](https://github.com/NCAR/summa/releases);
    * If you want the latest and greatest (and potentially erroneous), download a copy of the [development branch](https://github.com/ncar/summa/tree/develop) (or clone it);
    * If you may want to do SUMMA development, then fork the repo on github and start editing your own copy.

    Note that you will not be able to contribute to the main SUMMA repo directly. If you are seriously interested in contributing, spend a little time learning git. It will be useful anyway. For more information about working with the SUMMA code, please see the following documents:

    * [SUMMA and Git](../development/SUMMA_and_git.md)
    * [Git Workflow for SUMMA](../development/SUMMA_git_workflow.md)
    * [SUMMA Coding Conventions](../development/SUMMA_coding_conventions.md)

Once you have all the above, you can compile SUMMA using the following steps:

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

Continue reading [SUMMA configuration](../configuration/SUMMA_configuration.md) to learn more about how to configure SUMMA for your application. We strongly recommend that you get the [test applications](SUMMA_test_cases.md) to help you get started.
