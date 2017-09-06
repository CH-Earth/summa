# Instructions for installing SUMMA on Mac OS X

> The following has been tested on a MacBook Pro and iMac, both using OS X v. 10.12.6 (Sierra) and with [MacPorts](http://www.macports.org) as the OS X package manager. It should also work on older versions of OS X or with other package managers such as [fink](http://www.finkproject.org) or [homebrew](http://brew.sh), but you will undoubtedly need to make some modifications.

This document is mainly for people who want to use SUMMA in their modeling project, rather than contribute or change the SUMMA source code. If you want to contribute to SUMMA, , please see the following documents:

* [SUMMA and Git](../development/SUMMA_and_git.md)
* [Git Workflow for SUMMA](../development/SUMMA_git_workflow.md)
* [SUMMA Coding Conventions](../development/SUMMA_coding_conventions.md)


In the following I will assume that you don't have a Fortran compiler or NetCDF installed. If you do, just scroll down.

 1. Install [MacPorts](http://www.macports.org). This is a package manager for OS X (Darwin), which will let you install many software packages that are typically used in a Unix-like environment (which is what your Mac provides). MacPorts has detailed instructions on how it should be installed for various versions of OS X. Just follow the instructions (you will be asked to install Xcode, etcetera). If you know nothing about Unix/Linux/Darwin, now is the time to learn.

 1. Open a terminal (Applications --> Utilities --> Terminal.app)

 1. Use MacPorts to install NetCDF by typing the following on the command-line

    `sudo port install netcdf`

    You will be asked to provide your password. Note that you can only `sudo` (execute a command as a superuser) if you have administrative privileges for your machine (i.e., you can install software, etc.). If you do not have administrative privileges, talk to the person who maintains your computer. A whole bunch of other packages (dependencies) will be installed as well. Just let MacPorts do its work, this is why you are using a package manager. It may take a while the first time you do this. Make sure that the NetCDF is compiled using `gcc6` or higher. You can specify different variants of the package to be installed. See the MacPorts documentation for details.

 1. Install the Fortran version of the NetCDF library

    `sudo port install netcdf-fortran`

    Note what version of `gcc` was used to compile the NetCDF Fortran library. If you did not see it, type

    `sudo port info netcdf-fortran`

    The gcc compiler will be listed under `Build Dependencies`. For example, `gcc6`, which is what I'll assume for now. However, make sure to use the correct version of gcc (i.e. the one used to compile the netcdf-fortran library) when you compile SUMMA.

 1. Install the correct version of gcc (here we'll assume gcc6)

    `sudo port install gcc6`

 1. Note that MacPorts typically installs everything in the `/opt/local` directory. Make sure you now have NetCDF and gfortran installed:

    `ls /opt/local/lib/*netcdf*`

    You should have one or more `libnetcdf.*` (C-version of NetCDF library) and `libnetcdff.*` (Fortran version of NetCDF library) files

    `ls /opt/local/bin/*fortran*`

    You should have one or more `gfortran-*` files. If you installed `gcc6`, the gfortran executable will be `/opt/local/bin/gfortran-mp-6`. Since MacPorts will have modified your path so that `/opt/local/bin` is part of that, you should be able to invoke the compiler by typing

    `gfortran-mp-6`

    and the results should be

    `gfortran-mp-6: fatal error: no input files`
    `compilation terminated.`

    Note that you may also have a symbolic link named `/opt/local/bin/gfortran`. Make sure that that link points to the correct executable (e.g. `/opt/local/bin/gfortran-mp-6`). If so, then you should be able to invoke the correct version of the fortran compile by simply typing `gfortran`.

 1. While you are at it, there are a number of other packages that would be useful to install, in particular

    * `sudo port install ncview` : to visualize NetCDF files
    * `sudo port install nco`    : to manipulate NetCDF files, see the [NCO homepage](http://nco.sourceforge.net)
    * `sudo port install cdo`    : to manipulate NetCDF files, see the [CDO homepage](https://code.mpimet.mpg.de/projects/cdo/)

 1. Now obtain the SUMMA source code from the [SUMMA source code repository](https://github.com/NCAR/summa). You may just want to download the latest tagged release. Unless you are planning to contribute to the source code, there is no need to clone or fork the repository.

 1. Untar or unzip the archive, then go to the `summa/build` directory and follow the instructions in the [SUMMA installation](SUMMA_installation.md) page. If you are using MacPorts, the `FC_ENV` can be set to `gfortran-6-macports`.
