#!/bin/bash
  
# Build nextgen on Mac, from ngen directory put this one directory up and run this as ../build_ngen.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# activate correct python environment, here is an example with conda environment named ngen
# NOTE: ngen does not support numpy>=2.0, so this env must have numpy<2 (e.g. numpy 1.26.x)
: "${PYNGEN_CONDA_ENV:=ngen}"
# try common conda install locations; adjust if your conda is elsewhere
if [ -f "${HOME}/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/opt/anaconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)" || true
fi
# activate env if possible (non-fatal)
if command -v conda >/dev/null 2>&1; then
  conda activate "${PYNGEN_CONDA_ENV}" || true
fi
# fallback: allow overriding python executable explicitly
: "${NGEN_PYTHON_EXECUTABLE:=$(which python 2>/dev/null || echo /usr/bin/python3)}"
# root of the active python environment (asked of the interpreter itself so a stale
# VIRTUAL_ENV/CONDA_PREFIX can't mislead it); used as a hint for find_package(Python)
: "${NGEN_PYTHON_ROOT:=$("${NGEN_PYTHON_EXECUTABLE}" -c 'import sys; print(sys.prefix)' 2>/dev/null || dirname "$(dirname "${NGEN_PYTHON_EXECUTABLE}")")}"
# ngen does not support numpy>=2.0; verify the active env has numpy<2
"${NGEN_PYTHON_EXECUTABLE}" - <<'PY' || { echo "ERROR: need numpy<2 in the '${PYNGEN_CONDA_ENV}' env (e.g. conda install 'numpy<2')"; exit 1; }
import sys
from packaging.version import Version
import numpy as np
sys.exit(0 if Version(np.__version__) < Version("2.0") else 1)
PY
"${NGEN_PYTHON_EXECUTABLE}" -c 'import numpy as np; print("Using NumPy:", np.__version__)'

# Mac Example using MacPorts:
export CC=/opt/local/bin/gcc
export CXX=/opt/local/bin/g++
export FC=/opt/local/bin/gfortran
# C/C++ compiler to use for the ngen framework itself (set below, before the ngen
# cmake call).  MacPorts/Homebrew GCC cannot compile against the macOS >=15 / 26 SDK
# system headers from C++ (the <mach/*> headers fail with
# "expected constructor, destructor, or type conversion before '(' token"), so ngen
# and its C/C++ extern modules are built with Apple clang.  Fortran stays on gfortran.
# SUMMA and iso_c_fortran_bmi are plain Fortran/C and still build with GCC above; they
# connect to ngen only through the C-ABI BMI interface, so the mix is safe.
: "${NGEN_CC:=/usr/bin/clang}"
: "${NGEN_CXX:=/usr/bin/clang++}"

#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use
export C_INCLUDE_PATH=/opt/local/include
export CPLUS_INCLUDE_PATH=/opt/local/include
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=../../../sundials/instdir/                # will not be used if -DUSE_SUNDIALS=OFF

# Build SUMMA NGEN below, may wish to turn -DUSE_SUNDIALS=ON (must install Sundials first)
cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DUSE_NEXTGEN=ON -DUSE_SUNDIALS=OFF -DSPECIFY_LAPACK_LINKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build extern/summa/cmake_build --target all -j

# Build the ngen framework and its C/C++ extern modules with Apple clang (see note above).
export CC="${NGEN_CC}"
export CXX="${NGEN_CXX}"
cmake -S . -B cmake_build -DBoost_INCLUDE_DIR=/opt/local/libexec/boost/1.81/include \
    -DCMAKE_C_COMPILER="${NGEN_CC}"              \
    -DCMAKE_CXX_COMPILER="${NGEN_CXX}"           \
    -DCMAKE_Fortran_COMPILER="${FC}"             \
    -DPython_EXECUTABLE="${NGEN_PYTHON_EXECUTABLE}" \
    -DPython_ROOT_DIR="${NGEN_PYTHON_ROOT}"      \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo            \
    -DNGEN_IS_MAIN_PROJECT=ON                    \
    -DNGEN_WITH_MPI:BOOL=OFF                     \
    -DNGEN_WITH_NETCDF:BOOL=ON                   \
    -DNGEN_WITH_SQLITE:BOOL=ON                   \
    -DNGEN_WITH_UDUNITS:BOOL=ON                  \
    -DNGEN_WITH_BMI_FORTRAN:BOOL=ON              \
    -DNGEN_WITH_BMI_C:BOOL=ON                    \
    -DNGEN_WITH_PYTHON:BOOL=ON                   \
    -DNGEN_WITH_ROUTING:BOOL=ON                  \
    -DNGEN_WITH_TESTS:BOOL=ON                    \
    -DNGEN_QUIET:BOOL=ON                         \
    -DNGEN_WITH_EXTERN_ALL:BOOL=ON
    
# Comments on above choices, and defaults
#    -DCMAKE_BUILD_TYPE=RelWithDebInfo:  to be able to run in gdb change to -DCMAKE_BUILD_TYPE=Debug
#    -DNGEN_IS_MAIN_PROJECT=ON        :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_MPI:BOOL=OFF         :  may want to turn this ON as well as uncommenting "make -j 8 -C cmake_build"
#    -DNGEN_WITH_NETCDF:BOOL=ON       :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_SQLITE:BOOL=ON       :  must be BOOL=ON if planning to use GeoPackages (and not just geojsons)
#    -DNGEN_WITH_UDUNITS:BOOL=ON      :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_BMI_FORTRAN:BOOL=ON  :  must be BOOL=ON to build SUMMA NGEN
#    -DNGEN_WITH_BMI_C:BOOL=ON        :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_PYTHON:BOOL=ON       :  must be BOOL=ON for DNGEN_WITH_EXTERN_ALL:BOOL=ON
#    -DNGEN_WITH_ROUTING:BOOL=ON      :  must have DNGEN_WITH_PYTHON:BOOL=ON for this to be ON
#    -DNGEN_WITH_TESTS:BOOL=ON        :  must have DNGEN_WITH_EXTERN_ALL:BOOL=ON for this to be ON
#    -DNGEN_QUIET:BOOL=ON             :  may want turn to this OFF, especially if debugging
#    -DNGEN_WITH_EXTERN_ALL:BOOL=ON   :  these submodules are not used with SUMMA, you may turn this off you don't want to use them

# make -j 8 -C cmake_build    # build w/ 8 parallel jobs, if uncomment then comment the next line and use DNGEN_WITH_MPI:BOOL=ON
make -C cmake_build
