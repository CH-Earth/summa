#!/bin/bash
  
# Build nextgen on HPC, from ngen directory put this one directory up and run this as ../build_ngen.cluster.bash
# Load modules, example on Anvil
module load r/4.4.1
module load gcc/14.2.0
module load openmpi/4.1.6
module load gdal/3.10.0
module load conda/2024.09
module load openblas/0.3.17
module load netcdf-fortran/4.5.3
module load udunits/2.2.28
module load boost/1.86.0
module load sqlite/3.46.0-ayg27dg

# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# activate correct python environment, here is an example with conda environment named venv installed from SYMFLUENCE
: "${NGEN_CONDA_ENV:=venv}"
source ${HOME}/Symfluence/SYMFLUENCE/${NGEN_CONDA_ENV}/bin/activate
# fallback: allow overriding python executable explicitly
: "${NGEN_PYTHON_EXECUTABLE:=$(which python 2>/dev/null || echo /usr/bin/python3)}"
# root of the active python environment (asked of the interpreter itself so a stale
# VIRTUAL_ENV/CONDA_PREFIX can't mislead it); used as a hint for find_package(Python)
: "${NGEN_PYTHON_ROOT:=$("${NGEN_PYTHON_EXECUTABLE}" -c 'import sys; print(sys.prefix)' 2>/dev/null || dirname "$(dirname "${NGEN_PYTHON_EXECUTABLE}")")}"
export PYTHONNOUSERSITE=1
python -m pip install --upgrade "pip<24.1" >/dev/null 2>&1 || true
# Next Gen does not support numpy>=2.0; make sure the active env has numpy<2
python - <<'PY' || (python -m pip install "numpy<2" "setuptools<70" && true)
from packaging.version import Version
import numpy as np
assert Version(np.__version__) < Version("2.0")
PY
python - <<'PY'
import numpy as np
print("Using NumPy:", np.__version__)
PY

#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use

export SUNDIALS_DIR="$CMAKE_PREFIX_PATH:$HOME/Summa-Actors/utils/dependencies/install/sundials/"

# Build SUMMA NGEN below
cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
cmake --build extern/iso_c_fortran_bmi/cmake_build --target all

cmake -B extern/summa/cmake_build -S extern/summa -DUSE_NEXTGEN=ON -DUSE_SUNDIALS=OFF -DSPECIFY_LAPACK_LINKS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build extern/summa/cmake_build --target all -j

cmake -S . -B cmake_build \
    -DPython_EXECUTABLE="${NGEN_PYTHON_EXECUTABLE}" \
    -DPython_ROOT_DIR="${NGEN_PYTHON_ROOT}"      \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo            \
    -DNGEN_IS_MAIN_PROJECT=ON                    \
    -DNGEN_WITH_MPI:BOOL=ON                      \
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

