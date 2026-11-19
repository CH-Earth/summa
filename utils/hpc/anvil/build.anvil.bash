#!/usr/bin/env bash
set -euo pipefail

# Compile SUMMA for Anvil's GCC/OpenMPI/NetCDF stack.
REPO="$HOME/models/summa"
source "$REPO/utils/hpc/anvil/env_summa.anvil.bash"

if [[ ! -d "$SUNDIALS_DIR" ]]; then
  echo "ERROR: SUNDIALS is not installed at $SUNDIALS_DIR" >&2
  echo 'Build/install a compatible SUNDIALS release on Anvil first (see README).' >&2
  exit 1
fi
command -v mpifort >/dev/null || { echo 'ERROR: mpifort missing' >&2; exit 1; }
command -v nf-config >/dev/null || { echo 'ERROR: nf-config missing' >&2; exit 1; }
command -v cmake >/dev/null || { echo 'ERROR: cmake missing' >&2; exit 1; }

printf 'MPI Fortran: %s\n' "$(command -v mpifort)"
printf 'NetCDF Fortran: %s\n' "$(nf-config --version)"
printf 'SUNDIALS_DIR: %s\n' "$SUNDIALS_DIR"
module list

# Make the local SUNDIALS install discoverable without replacing other CMake paths.
export CMAKE_PREFIX_PATH="$SUNDIALS_DIR${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"

cd "$REPO/build/cmake"
rm -rf ../cmake_build
cmake -B ../cmake_build -S ../. \
  -DUSE_MPI=ON \
  -DUSE_SUNDIALS=ON \
  -DUSE_MIZUROUTE=ON \
  -DSPECIFY_LAPACK_LINKS=OFF \
  -DCMAKE_BUILD_TYPE=Release

cmake --build ../cmake_build --target all -j "${SUMMA_BUILD_JOBS:-8}"

test -x "$REPO/bin/summa_sundials_mizuroute_opt.exe" || {
  echo 'ERROR: Expected calibration executable was not created.' >&2
  exit 1
}
echo 'Anvil SUMMA build complete.'
