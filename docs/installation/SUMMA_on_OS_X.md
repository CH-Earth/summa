# Instructions for installing SUMMA on Mac OS X

Two slightly different options have been tested to work. The first uses MacPorts as package manager. The second uses Homebrew as package manager.

This document is mainly for people who want to use SUMMA in their modeling project, rather than contribute or change the SUMMA source code. If you want to contribute to SUMMA, please see the following documents:

* [SUMMA and Git](../development/SUMMA_and_git.md)
* [Git Workflow for SUMMA](../development/SUMMA_git_workflow.md)
* [SUMMA Coding Conventions](../development/SUMMA_coding_conventions.md)

## Instructions with MacPorts
> The following has been tested on a MacBook Pro and iMac, both using OS X v. 10.12.6 (Sierra) and with [MacPorts](http://www.macports.org) as the OS X package manager. It should also work on older versions of OS X or with other package managers such as [fink](http://www.finkproject.org) or [homebrew](http://brew.sh), but you will undoubtedly need to make some modifications.

In the following I will assume that you don't have a Fortran compiler or NetCDF installed. If you do, just scroll down.

### Step 1: install MacPorts
 1. Install [MacPorts](http://www.macports.org). This is a package manager for OS X (Darwin), which will let you install many software packages that are typically used in a Unix-like environment (which is what your Mac provides). MacPorts has detailed instructions on how it should be installed for various versions of OS X. Just follow the instructions (you will be asked to install Xcode, etcetera). If you know nothing about Unix/Linux/Darwin, now is the time to learn.

### Step 2: install packages

 1. Open a terminal (Applications --> Utilities --> Terminal.app)

 2. Use MacPorts to install NetCDF by typing the following on the command-line

    `sudo port install netcdf`

    You will be asked to provide your password. Note that you can only `sudo` (execute a command as a superuser) if you have administrative privileges for your machine (i.e., you can install software, etc.). If you do not have administrative privileges, talk to the person who maintains your computer. A whole bunch of other packages (dependencies) will be installed as well. Just let MacPorts do its work, this is why you are using a package manager. It may take a while the first time you do this. Make sure that the NetCDF is compiled using `gcc6` or higher. You can specify different variants of the package to be installed. See the MacPorts documentation for details.

 3. Install the Fortran version of the NetCDF library

    `sudo port install netcdf-fortran`

    Note what version of `gcc` was used to compile the NetCDF Fortran library. If you did not see it, type

    `sudo port info netcdf-fortran`

    The gcc compiler will be listed under `Build Dependencies`. For example, `gcc6`, which is what I'll assume for now. However, make sure to use the correct version of gcc (i.e. the one used to compile the netcdf-fortran library) when you compile SUMMA.

 4. Install the correct version of gcc (here we'll assume gcc6)

    `sudo port install gcc6`

 5. Note that MacPorts typically installs everything in the `/opt/local` directory. Make sure you now have NetCDF and gfortran installed:

    `ls /opt/local/lib/*netcdf*`

    You should have one or more `libnetcdf.*` (C-version of NetCDF library) and `libnetcdff.*` (Fortran version of NetCDF library) files

    `ls /opt/local/bin/*fortran*`

    You should have one or more `gfortran-*` files. If you installed `gcc6`, the gfortran executable will be `/opt/local/bin/gfortran-mp-6`. Since MacPorts will have modified your path so that `/opt/local/bin` is part of that, you should be able to invoke the compiler by typing

    `gfortran-mp-6`

    and the results should be

    `gfortran-mp-6: fatal error: no input files`
    `compilation terminated.`

    Note that you may also have a symbolic link named `/opt/local/bin/gfortran`. Make sure that that link points to the correct executable (e.g. `/opt/local/bin/gfortran-mp-6`). If so, then you should be able to invoke the correct version of the fortran compile by simply typing `gfortran`.

 6. While you are at it, there are a number of other packages that would be useful to install, in particular

    * `sudo port install ncview` : to visualize NetCDF files
    * `sudo port install nco`    : to manipulate NetCDF files, see the [NCO homepage](http://nco.sourceforge.net)
    * `sudo port install cdo`    : to manipulate NetCDF files, see the [CDO homepage](https://code.mpimet.mpg.de/projects/cdo/)


### Step 3: Compile SUMMA
 1. Now obtain the SUMMA source code from the [SUMMA source code repository](https://github.com/NCAR/summa). You may just want to download the latest tagged release. Unless you are planning to contribute to the source code, there is no need to clone or fork the repository.

 2. Untar or unzip the archive, then go to the `summa/build` directory and follow the instructions in the [SUMMA installation](SUMMA_installation.md) page. If you are using MacPorts, the `FC_ENV` can be set to `gfortran-6-macports`.

## Instructions with Homebrew

 These instructions deviate from the standard SUMMA on OS X installation instructions in that Homebrew is used as an alternative
 to MacPorts as a package installer. This is largely due to Catalina changes to a users /opt/ folder.

 These instructions assume that you've already cloned the SUMMA repository.
 If you haven't, follow the initial [SUMMA Installation Instructions](https://summa.readthedocs.io/en/latest/installation/SUMMA_installation/)

### Step 1: Installation of Homebrew

 1. If you are already using MacPorts, it is recommended to uninstall it to avoid package issues.
 Uninstallation instructions are available [here](https://guide.macports.org/chunked/installing.macports.uninstalling.html).
 If you'd like to save the currently installed packages for MacPorts to file, use the command: `port installed > ports_installed.txt`

 2. Homebrew can then be installed following the instructions [here](https://brew.sh/)

### Step 2: Installation of Required packages

 1. Use the command as follows to install gcc and gfortran compilers.
 `brew install gcc`

 2. Use the command as follows to install netCDF.
 `brew install netcdf`

 Note that as a bonus, NetCDF Fortran is installed with this command and does not require action.
 It is installed as a [resource](https://github.com/Homebrew/homebrew-core/blob/HEAD/Formula/netcdf.rb), by default.

 3. Use the following commands to install lapack and openblas using Homebrew
 `brew install lapack`
 `brew install openblas`

### Step 3: Compile SUMMA
  1. Now obtain the SUMMA source code from the [SUMMA source code repository](https://github.com/NCAR/summa). You may just want to download the latest tagged release. Unless you are planning to contribute to the source code, there is no need to clone or fork the repository.

 2. Untar or unzip the archive, then go to the `summa/build` directory and follow the instructions in the [SUMMA installation](SUMMA_installation.md) page. If you are using Homebrew, you can also follow the following steps:

 3. Navigate to your local copy of the SUMMA directory and go to the build subdirectory;
 4. Make a copy of the Makefile, naming it Makefile.local for your own use.
 5. Update the Environment Variables in the Makefile.local following the instructions. What worked for me is:

 F_MASTER=$PATH_TO_SUMMA_GIT_CLONE$
 FC=gfortran
 FC_EXE=/usr/local/bin/gfortran
 INCLUDES=-I/usr/local/include -I/usr/local/opt/lapack/include -I/usr/local/opt/openblas/include
 LIBRARIES=-L/usr/local/lib -lnetcdff -L/usr/local/opt/lapack/lib -lblas -L/usr/local/opt/openblas/lib -lopenblas

 While Homebrew uses the /usr/local/ dir and this should be consistent, it is recommend to check these folders to verify the files exist.

 Navigating to the SUMMA base directory (i.e. F_MASTER, the root directory of the build folder), run the command

 6. `make -f build/Makefile.local`

 If the code compiles successfully, then the last line of output from the make process will tell you where the SUMMA executable is installed (it goes into summa/bin). Run summa.exe in that directory (you may need to provide the full path).
