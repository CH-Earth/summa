#!/bin/sh
# Fetch OpenWQ and build its local third-party dependencies (Armadillo, PhreeqcRM)
# into build/source/openwq/_deps/. Idempotent: re-running skips work already done.
#
# See README.md for the full picture and the SUMMA cmake configure line.
#
# Overridable via environment:
#   CXX, CC     C++/C compilers (MUST match the compiler SUMMA is built with)
#   DEPS_DIR    where to install the built dependencies   (default: ./_deps)
#   ARMA_VER    Armadillo version to download             (default: 14.2.3)
#   JOBS        parallel build jobs                       (default: 4)
#   H5_API      HDF5 legacy-API define for Armadillo      (default: H5_USE_110_API)
#   H5_PREFIX   prefix containing hdf5.h / libhdf5        (default: /opt/local)
set -eu

HERE=$(cd "$(dirname "$0")" && pwd)
CXX=${CXX:-/opt/local/bin/g++}
CC=${CC:-/opt/local/bin/gcc}
DEPS_DIR=${DEPS_DIR:-"$HERE/_deps"}
ARMA_VER=${ARMA_VER:-14.2.3}
JOBS=${JOBS:-4}
H5_API=${H5_API:-H5_USE_110_API}
H5_PREFIX=${H5_PREFIX:-/opt/local}

echo ">> compilers : CXX=$CXX  CC=$CC"
echo ">> deps dir  : $DEPS_DIR"

# ---------------------------------------------------------------------------
# OpenWQ C++ source
# ---------------------------------------------------------------------------
if [ ! -d "$HERE/openwq/src" ]; then
    echo ">> cloning OpenWQ (develop) into openwq/"
    git clone -b develop https://github.com/ue-hydro/openwq.git "$HERE/openwq"
else
    echo ">> openwq/ already present - skipping clone"
fi

mkdir -p "$DEPS_DIR"

# ---------------------------------------------------------------------------
# Armadillo  (built WITH ARMA_USE_HDF5 - OpenWQ needs arma::hdf5_name / hid_t)
# ---------------------------------------------------------------------------
ARMA_PREFIX="$DEPS_DIR/armadillo-install"
if [ ! -f "$ARMA_PREFIX/lib/libarmadillo.dylib" ] && [ ! -f "$ARMA_PREFIX/lib/libarmadillo.so" ]; then
    echo ">> building Armadillo $ARMA_VER"
    cd "$DEPS_DIR"
    if [ ! -d "armadillo-$ARMA_VER" ]; then
        curl -L -o "armadillo-$ARMA_VER.tar.xz" \
          "https://sourceforge.net/projects/arma/files/armadillo-$ARMA_VER.tar.xz/download"
        tar -xf "armadillo-$ARMA_VER.tar.xz"
    fi
    cd "armadillo-$ARMA_VER"
    # enable ARMA_USE_HDF5 in the config template (regenerated at cmake time)
    perl -0pi -e 's{^// #define ARMA_USE_HDF5$}{#define ARMA_USE_HDF5}m' \
        include/armadillo_bits/config.hpp.cmake
    rm -rf build
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_PREFIX_PATH="$H5_PREFIX" \
        -DCMAKE_CXX_FLAGS="-D$H5_API -I$H5_PREFIX/include" \
        -DCMAKE_SHARED_LINKER_FLAGS="-L$H5_PREFIX/lib -lhdf5" \
        -DCMAKE_INSTALL_PREFIX="$ARMA_PREFIX"
    cmake --build build -j"$JOBS"
    cmake --install build
else
    echo ">> Armadillo already installed at $ARMA_PREFIX - skipping"
fi

# ---------------------------------------------------------------------------
# PhreeqcRM
# ---------------------------------------------------------------------------
PHREEQC_PREFIX="$DEPS_DIR/phreeqcrm-install"
if [ ! -f "$PHREEQC_PREFIX/include/PhreeqcRM.h" ]; then
    echo ">> building PhreeqcRM"
    cd "$DEPS_DIR"
    [ -d phreeqcrm ] || git clone --depth 1 https://github.com/usgs-coupled/phreeqcrm.git phreeqcrm
    cd phreeqcrm
    rm -rf build
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" \
        -DBUILD_SHARED_LIBS=ON \
        -DPHREEQCRM_BUILD_TESTS=OFF -DPHREEQCRM_FORTRAN_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX="$PHREEQC_PREFIX"
    cmake --build build -j"$JOBS"
    cmake --install build
else
    echo ">> PhreeqcRM already installed at $PHREEQC_PREFIX - skipping"
fi

# ---------------------------------------------------------------------------
cat <<EOF

>> dependencies ready.

Now configure SUMMA from the build/ directory, e.g.:

  DEPS=$DEPS_DIR
  SUNDIALS=<path to your SUNDIALS install>

  FC=gfortran LIBRARY_LINKS="-framework Accelerate" \\
  cmake -B cmake_build_openwq -S . \\
    -DCMAKE_BUILD_TYPE=Release -DUSE_SUNDIALS=ON -DUSE_OPENWQ=ON -DSPECIFY_LAPACK_LINKS=ON \\
    -DSUNDIALS_DIR="\$SUNDIALS/lib/cmake/sundials" \\
    -DPhreeqcRM_DIR="\$DEPS/phreeqcrm-install/lib/cmake/PhreeqcRM" \\
    -DCMAKE_PREFIX_PATH="$H5_PREFIX;\$DEPS/armadillo-install;\$DEPS/phreeqcrm-install;\$SUNDIALS"

  cmake --build cmake_build_openwq -j$JOBS
EOF
