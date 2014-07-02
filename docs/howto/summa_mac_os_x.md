# Instructions for installing SUMMA on Mac OS X

> The following has been tested on a MacBook Pro and iMac, both using OS X v. 10.9.3 (Mavericks) and with [MacPorts](http://www.macports.org) as the OS X package manager. It should also work on older versions of OS X or with other package managers such as [fink](http://www.finkproject.org) or [homebrew](http://brew.sh), but you will undoubtedly need to make some modifications.

This document is mainly for people who want to use SUMMA in their modeling project, rather than contribute or change the SUMMA source code. If you want to contribute to SUMMA, take a look at the git related documents in the [docs/howto](docs/howto) directory.

In the following I will assume that you don't have a Fortran compiler or NetCDF installed. If you do, just scroll down.

 1. Install [MacPorts](http://www.macports.org). This is a package manager for OS X (Darwin), which will let you install many software packages that are typically used in a Unix-like environment (which is what your Mac provides). MacPorts has detailed instructions on how it should be installed for various versions of OS X. Just follow the instructions (you will be asked to install Xcode, etcetera). If you know nothing about Unix/Linux/Darwin, now is the time to learn.

 1. Open a terminal (Applications --> Utilities --> Terminal.app)

 1. Use MacPorts to install NetCDF

        sudo port install netcdf

    You will be asked to provide your password. Note that you can only `sudo` (execute a command as a superuser) if you have administrative privileges for your machine (i.e., you can install software, etc.). A whole bunch of other packages (dependencies) will be installed as well. Just let MacPorts do its work, this is why you are using a package manager. It may take a while the first time you do this.

 1. Install the Fortran version of the NetCDF library

        sudo port install netcdf-fortran

    Note what version of `gcc` was used to compile the NetCDF Fortran library. If you did not see it, type

        sudo port info netcdf-fortran

    The gcc compiler will be listed under 'Dependencies'. For example, `gcc48`, which is what we'll assume for now. However, make sure to use the correct version of gcc (i.e. the one used to compile the netcdf-fortran library)

 1. Install the correct version of gcc (here we'll assume gcc48)

        sudo port install gcc48

 1. Note that MacPorts typically installs everything in the `/opt/local` directory. Make sure you now have NetCDF and gfortran installed:

        ls /opt/local/lib/*netcdf*

    You should have one or more `libnetcdf.*` (C-version of NetCDF library) and `libnetcdff.*` (Fortran version of NetCDF library) files

        ls /opt/local/bin/*fortran*

    You should have one or more `gfortran-*` files. If you installed `gcc48`, the gfortran executable will be `/opt/local/bin/gfortran-mp-4.8`. Since MacPorts will have modified your path so that `/opt/local/bin` is part of that, you should be able to invoke the compiler by typing

        gfortran-mp-4.8

    and the results should be

        gfortran-mp-4.8: fatal error: no input files
        compilation terminated.

 1. While you are at it, there are a number of other packages that would be useful to install, in particular

        sudo port install ncview

 1. Now obtain the SUMMA source code from the [SUMMA source code repository](https://github.com/UW-Hydro/summa). You may just want to download the latest tagged release. Since you are not planning to contribute to the source code, there is no need to clone or fork the repository.

 1. Untar or unzip the archive, then go to the `summa/build` directory.

 1. Modify the Makefile. At the very least you will need to modify `F_MASTER`, the root directory for summa and `FC`, the compiler. For example, if your `summa` directory is in `/home/user/summa` and the compiler is `gfortran-mp-4.8`, then simply set the following:

        F_MASTER = /home/user/summa
        FC = gfortran-mp-4.8

    Makefiles can be finicky. Make sure there are no trailing spaces on these lines. You want to do the editing in a text editor (vi, vim, emacs, sublime, nano, textedit, bbedit, whatever you prefer, just not something like MS Word or Pages). Note that `F_MASTER` is specified near the top of the file and `FC` about halfway down

 1. Check that all variables in the Makefile are set correctly

        make check

    Inspect the variables and make sure that they make sense. If not, modify the Makefile further. The `NCDF_PATH` on OS X is `/opt/local` when you use MacPorts

 1. If everything looks OK, then

        make

    This will build `summa.exe` and move it to the `summa/bin` directory if it is successful

 1. Pay attention to the `make` output. You may need to set some environment variables (`LD_LIBRARY_PATH` in particular) to support dynamic linking

 1. Run `summa.exe`. If all goes well, you should get an error message that looks something like:

    ```
    231714.555
    1st command-line argument missing, expect text string defining the output file suffix
    ```

If you get this far then SUMMA is installed correctly and functional. See the [SUMMA web site](http://www.ral.ucar.edu/projects/summa) at NCAR for sample data sets and additional information.
