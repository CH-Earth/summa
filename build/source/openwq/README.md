# OpenWQ Integration

This directory contains the source code that enables Summa to couple to openWQ.

To compile with openWQ support, you need to do the following steps:
  - clone the openWQ repository: `git clone -b develop https://github.com/ue-hydro/openwq.git`
    - The above needs to be done in build/source/openwq
  - To compile Summa-OpenWQ, compile summa normally with CMake, but add the flag `-DENABLE_OPENWQ=ON`
    - Compiling with openWQ support works for both Sundials and non-sundials builds of SUMMA. 