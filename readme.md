#SUMMA

This is the source code repository for the **Structure for Unifying Multiple Modeling Alternatives** or **SUMMA**. More information about SUMMA, including publications, test data sets, and sample applications can be found on the [SUMMA web site](http://www.ral.ucar.edu/projects/summa) at NCAR.

##Description
SUMMA is a hydrologic modeling aproach that is built on a common set of governing equations and a common numerical solver, which together constitute the  “structural core” of the model. Different modeling approaches can then be implemented within the structural core, enabling a controlled and systematic analysis of alternative modeling options, and providing insight for future model development.

The important modeling features are:

 1. The formulation of the governing model equations is cleanly separated from their numerical solution;

 1. Different model representations of physical processes (in particular, different flux parameterizations) can be used within a common set of governing equations; and

 1. The physical processes can be organized in different spatial configurations, including model elements of different shape and connectivity (e.g., nested multi-scale grids and HRUs).

SUMMA can be used to configure a wide range of hydrological model alternatives. We anticipate that systematic model analysis will help researchers and practitioners understand reasons for inter-model differences in model behavior, and, when applied across a large sample of catchments, may provide insights on the dominance of different physical processes and regional variability in the suitability of different modeling approaches. An important application of SUMMA is selecting specific physics options to reproduce the behavior of existing models – these applications of *“model mimicry”* can be used to define reference (benchmark) cases in structured model comparison experiments, and can help diagnose weaknesses of individual models in different hydroclimatic regimes.

##Credits
SUMMA's initial implementation is described in two papers that have been submitted to [Water Resources Research](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1944-7973). If you use SUMMA, please credit these two publications.

 * Clark, M. P., B. Nijssen, J. Lundquist, D. Kavetski, D. Rupp, R. Woods, E. Gutmann, A. Wood, L. Brekke, J. Arnold, D. Gochis, R. Rasmussen, 2014: A unified approach to hydrologic modeling: Part 1. Model structure. *Water Resources Research*, in review.

 * Clark, M. P., B. Nijssen, J. Lundquist, D. Kavetski, D. Rupp, R. Woods, E. Gutmann, A. Wood, D. Gochis, R. Rasmussen, D. Tarboton, V. Mahat, G. Flerchinger, D. Marks 2014: A unified approach to hydrologic modeling: Part 2. Comparison of alternative process representations. *Water Resources Research*, in review.

##Installation

We have successfully installed SUMMA on a number of Unix-like (*nix) operating systems, including Linux and Darwin (Mac OS X). To compile SUMMA, you will need:

 * a Fortran compiler. We have successfully used the intel Fortran compiler (`ifort`), the Portland Group Fortran 90 compiler (`pgf90`), and the GNU Fortran compiler (`gfortran`), the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile SUMMA with other Fortran compilers as well. If you do, please let us know and send us a copy of your Makefile.

    If you do not have a Fortran compiler, ou can install `gfortran` for free. The easiest way is to use a apackage manager. Which package manager depends on your *nix flavor. On OS X, you can use any of the free OS X package managers, including [MacPorts](http://www.macports.org), [fink](http://www.finkproject.org), or [homebrew](http://brew.sh). Note that `gfortran` is installed as part of the `gcc` compiler suite.

 * the NetCDF libraries. [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. They are widely used in the hydrometeorological community and eventually all SUMMA I/O will use NetCDF. Most *nix package managers (including  include a NetCDF port. Note that you need to ensure that:

    * You have NetCDF version 4.x;
    * The NetCDF libraries are compiled with the same compiler as you plan to use for compiling SUMMA; and
    * You have the NetCDF Fortran library installed (`libnetcdff.*`) and not just the C-version.


 * a copy of the SUMMA source code from [this repo](https://github.com/UW-Hydro/summa). You have a number of options:

    * If you just want to use the latest stable release of SUMMA, then simply look for the most recent tag;
    * If you want the latest and greatest (and potentially erroneous), download a copy of the master branch (or clone it);
    * If you may want to do SUMMA development, then fork the repo on github and start editing your own copy.

    Note that you will not be able to contribute to the main SUMMA repo directly. If you are seriously interested in contributing, spend a little time learning git. It will be useful anyway. For more information about working with the SUMMA code, please see the following documents:

    * [SUMMA and Git](https://github.com/UW-Hydro/summa/docs/summa_and_git_howto.md)
    * [Git workflow for SUMMA](https://github.com/UW-Hydro/summa/docs/summa_git_workflow.md)
    * [Git commands](https://github.com/UW-Hydro/summa/docs/git_howto.md)

Once you have all the above, you can compile SUMMA using the following steps:

 1. Navigate to your local copy of the SUMMA directory and go to the `build` subdirectory;

 1. Edit the `Makefile`. At the very least, you will need to set `F_MASTER` and `FC`. You may also need to set `NCDF_PATH` and you may need to add some extra entries if you are using a different Fortran compiler or your setup is different;

 1. Type `make`. If all goes well, this will build SUMMA and move the executable `summa.exe` to the `bin` directory;

 1. Pay attention to the `make` output. You may need to set some environment variables (`LD_LIBRARY_PATH` in particular) to support dynamic linking;

 1. Run `summa.exe`. If all goes well, you should get an error message that looks something like:

    ```
    231714.555
    1st command-line argument missing, expect text string defining the output file suffix
    ```

If you get this far then SUMMA is installed correctly and functional.

