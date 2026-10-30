#!/bin/bash
# Build SUMMA develop's own calibration program (the Fortran async DDS: summa_driver_opt.f90,
# summa_parameter_sampling.f90, parameter_search.f90) from the SAME develop checkout as libsumma.so,
# into a SEPARATE build folder. The existing build/cmake_build/libsumma.so is not touched.
#   result: summa-develop/bin/summa[_sundials]_mizuroute_opt.exe  (calibration, needs MPI)
#           summa-develop/bin/summa[_sundials]_mizuroute.exe      (one model run, used for the check)
# Run on an Anvil login node:  bash my_work_athabasca/scripts/fortran_dds/build_develop_opt.sh
set -euo pipefail
DEV=<WORK_ROOT>/summa-develop
module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3
cd "$DEV"
echo "develop at: $(git log -1 --format='%h %ad %s')"

# same SUNDIALS setting as the libsumma.so build, so the model physics match J()
SUN=OFF
grep -q "^USE_SUNDIALS:BOOL=ON" build/cmake_build/CMakeCache.txt 2>/dev/null && SUN=ON
SFX=""; [ "$SUN" = ON ] && SFX="_sundials"
echo "USE_SUNDIALS=$SUN (copied from build/cmake_build)"

cmake -B build/cmake_build_opt -S build \
    -DUSE_SUNDIALS=$SUN -DUSE_MPI=ON -DUSE_MIZUROUTE=ON -DUSE_NEXTGEN=OFF -DUSE_OPENWQ=OFF \
    -DSPECIFY_LAPACK_LINKS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake_build_opt --target summa${SFX}_mizuroute_opt.exe summa${SFX}_mizuroute.exe -j 8

ls -l bin/summa${SFX}_mizuroute_opt.exe bin/summa${SFX}_mizuroute.exe
echo "OK: Fortran DDS calibration program built"
