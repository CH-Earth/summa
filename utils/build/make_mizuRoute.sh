#!/usr/bin/env bash
set -e

# always run from the SUMMA root
SUMMA_ROOT="$(pwd)"

MIZU_ROOT="${SUMMA_ROOT}/external/mizuRoute"
MIZU_BUILD="${MIZU_ROOT}/route/build"
MIZU_BIN="${MIZU_ROOT}/route/bin"

PIO_DIR="${MIZU_ROOT}/libraries/parallelio"
PIO_TAG="pio2_6_6"

OUT_BIN="${SUMMA_ROOT}/bin"

mkdir -p "${OUT_BIN}"
mkdir -p "${MIZU_ROOT}/libraries"

# ----------------------------------------------------------------------
# ParallelIO dependency
# ----------------------------------------------------------------------

if [ ! -d "${PIO_DIR}/.git" ]; then
    git clone --branch "${PIO_TAG}" \
              https://github.com/NCAR/ParallelIO \
              "${PIO_DIR}"
fi

# ----------------------------------------------------------------------
# Clean previous mizuRoute build
# ----------------------------------------------------------------------

rm -rf "${MIZU_BUILD}/lib"
rm -f  "${MIZU_BUILD}/route_runoff"
rm -rf "${MIZU_BUILD}/route_runoff.dSYM"

# restore tracked build infrastructure removed above
git -C "${MIZU_ROOT}" restore ${MIZU_BUILD}/lib

# make sure the PIO build directory exists
mkdir -p "${MIZU_BUILD}/lib/piolib"

# ----------------------------------------------------------------------
# Build standalone mizuRoute
# ----------------------------------------------------------------------

cd "${MIZU_BUILD}"

export BLDDIR="$(pwd)/../"
export FC=gnu
export FC_EXE=mpif90
export CC=mpicc

export NCDF_PATH=/opt/homebrew/
export PNETCDF_PATH=/opt/homebrew/opt/pnetcdf

make \
    FC="${FC}" \
    FC_EXE="${FC_EXE}" \
    F_MASTER="${BLDDIR}" \
    NCDF_PATH="${NCDF_PATH}" \
    PNETCDF_PATH="${PNETCDF_PATH}" \
    EXE=route_runoff

# ----------------------------------------------------------------------
# Copy executable to SUMMA bin directory
# ----------------------------------------------------------------------

cd "${SUMMA_ROOT}"

cp "${MIZU_BIN}/route_runoff" "${OUT_BIN}/route_runoff"

# remove macOS debug-symbol bundle created inside the submodule
rm -rf "${MIZU_BUILD}/route_runoff.dSYM"

echo
echo "Built ${OUT_BIN}/route_runoff"
