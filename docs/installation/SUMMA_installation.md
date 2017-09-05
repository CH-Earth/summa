# SUMMA Installation

We have successfully installed SUMMA on a number of Unix-like (*nix) operating systems, including Linux and Darwin (Mac OS X). To compile SUMMA, you will need:

 * a Fortran compiler. We have successfully used the intel Fortran compiler (`ifort`) and the GNU Fortran compiler (`gfortran`, version 6 or higher), the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile SUMMA with other Fortran compilers as well.

    If you do not have a Fortran compiler, you can install `gfortran` for free. The easiest way is to use a package manager. Which package manager depends on your *nix flavor. On OS X, you can use any of the free OS X package managers, including [MacPorts](http://www.macports.org), [fink](http://www.finkproject.org), or [homebrew](http://brew.sh). Note that `gfortran` is installed as part of the `gcc` compiler suite.

 * the NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. They are widely used in the hydrometeorological community and eventually all SUMMA I/O will use NetCDF. Most *nix package managers include a NetCDF port. Note that you need to ensure that:

    * You have NetCDF version 4.x;
    * The NetCDF libraries are compiled with the same compiler as you plan to use for compiling SUMMA; and
    * You have the NetCDF Fortran library installed (`libnetcdff.*`) and not just the C-version.

 * the LAPACK — Linear Algebra PACKage library. [LAPACK](http://www.netlib.org/lapack/) provides a series of routines for linear algebra operations, including matrix solvers. How to install the library depends on your *nix variant and is not covered here. For example, on OS X you will get all the necessary LAPACK routines by installing the ATLAS software (again, this is easiest using a package manager; note that ATLAS can take many hours to build the first time when you install it).

 * The [ATLAS](http://math-atlas.sourceforge.net/) (Automatically Tuned Linear Algebra Software) library. Note that this is required on OS X using the gfortran compiler to be able to use LAPACK. It's not clear that this is used on other *nix machines.

 * a copy of the SUMMA source code from [this repo](https://github.com/NCAR/summa). You have a number of options:

    * If you just want to use the latest stable release of SUMMA, then simply look for the most recent tag;
    * If you want the latest and greatest (and potentially erroneous), download a copy of the master branch (or clone it);
    * If you may want to do SUMMA development, then fork the repo on github and start editing your own copy.

    Note that you will not be able to contribute to the main SUMMA repo directly. If you are seriously interested in contributing, spend a little time learning git. It will be useful anyway. For more information about working with the SUMMA code, please see the following documents:

    * [SUMMA and Git](https://github.com/NCAR/summa/blob/master/docs/howto/summa_and_git_howto.md)
    * [Git workflow for SUMMA](https://github.com/NCAR/summa/blob/master/docs/howto/summa_git_workflow.md)
    * [Git commands](https://github.com/NCAR/summa/blob/master/docs/howto/git_howto.md)
    * [SUMMA coding conventions](https://github.com/NCAR/summa/blob/master/docs/howto/summa_coding_conventions.md)

Once you have all the above, you can compile SUMMA using the following steps:

 1. Navigate to your local copy of the SUMMA directory and go to the `build` subdirectory;

 1. Edit part 0 of the `Makefile`. At the very least, you will need to set `F_MASTER` and `FC`. You may also need to set `NCDF_PATH` and `LAPK_PATH` and you may need to add some extra entries if you are using a different Fortran compiler or your setup is different (if someone wants to contribute an actual `configure` script that would be great);

 1. Type `make`. If all goes well, this will build SUMMA and move the executable `summa.exe` to the `bin` directory. You may get some warnings (depending on your compiler settings), but you should not get any errors;

 1. Pay attention to the `make` output. You may need to set some environment variables (`LD_LIBRARY_PATH` in particular) to support dynamic linking;

 1. Run `summa.exe`. If all goes well, you should get an error message that looks something like:

    ```
    Usage: summa.exe master_file [-s file_suffix] [-g startGRU countGRU] [-c checkHRU]
      summa.exe   -- summa executable
      file_suffix -- text string defining the output file suffix
      master_file -- path/name of master file
      startGRU    -- the index of the first GRU for a parallelization run
      countGRU    -- the number of GRUs for a parallelization run
      checkHRU    -- the index of the HRU for a single HRU run
    ```

If you get this far then SUMMA is installed correctly and functional.

Continue reading [SUMMA configuration](https://github.com/NCAR/summa/blob/master/docs/howto/summa_configuration.md) to learn more about how to configure SUMMA for your application. We strongly recommend that you get the [test applications](http://ral.ucar.edu/projects/summa) to help you get started.
