#!/bin/bash
# Build NCAR/summa develop with the SPOTPY C wrapper into a NEW folder, next to the
# existing summa-spotpy build (which is not touched).
#   result: <WORK_ROOT>/summa-develop/build/cmake_build/libsumma.so
# Run on an Anvil login node:  bash build_develop.sh
set -euo pipefail
ROOT=<WORK_ROOT>
OLD=$ROOT/summa-spotpy                 # current, validated build (reference)
NEW=$ROOT/summa-develop
PORT=$(cd "$(dirname "$0")" && pwd)     # folder holding summa_c_api.f90 (this script's folder)

module load gcc/14.2.0 openmpi/4.1.6 openblas/0.3.17 netcdf-fortran/4.5.3

# 1. get develop with its submodules (mizuRoute, toml-f)
if [ ! -d "$NEW/.git" ]; then
  git clone --branch develop --recurse-submodules https://github.com/NCAR/summa.git "$NEW"
fi
cd "$NEW"
echo "develop at: $(git log -1 --format='%h %ad %s')"

# 2. add the wrapper to SUMMA's shared library (libsumma.so)
cp "$PORT/summa_c_api.f90" build/source/driver/
grep -q summa_c_api.f90 build/source/CMakeLists.txt || \
  sed -i 's|    ${DRIVER_DIR}/summa_spinup.f90|    ${DRIVER_DIR}/summa_spinup.f90\n    ${DRIVER_DIR}/summa_c_api.f90|' build/source/CMakeLists.txt
grep -n "summa_c_api" build/source/CMakeLists.txt

# 3. use the same SUNDIALS setting as the current build, so the model physics match
SUN=OFF
if grep -q "^USE_SUNDIALS:BOOL=ON" "$OLD/build/cmake_build/CMakeCache.txt" 2>/dev/null; then SUN=ON; fi
echo "USE_SUNDIALS=$SUN (copied from $OLD)"

# 4. configure and build (Release, mizuRoute coupled; TOML reader is switched on automatically)
cmake -B build/cmake_build -S build \
    -DUSE_SUNDIALS=$SUN -DUSE_MPI=OFF -DUSE_NEXTGEN=OFF -DUSE_OPENWQ=OFF \
    -DUSE_MIZUROUTE=ON -DSPECIFY_LAPACK_LINKS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake_build --target summa -j 8

# 5. check the library and the exported symbol
LIB=$(find build/cmake_build -name "libsumma.so*" -type f | head -1)
ls -l "$LIB"
nm -D "$LIB" | grep -E " T summa_evaluate$" && echo "OK: summa_evaluate exported"
