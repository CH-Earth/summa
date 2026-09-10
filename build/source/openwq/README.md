# OpenWQ Integration

This directory contains the SUMMA-side coupling to [OpenWQ](https://github.com/ue-hydro/openwq):

| file | purpose |
|------|---------|
| `summa_openWQ.f90`            | passes SUMMA state volumes and inter-compartment fluxes to OpenWQ; walks the `gru -> hru -> dom -> var` spatial-domain data structures (every HRU is subdivided into upland / glacier / wetland domains, each an independent OpenWQ spatial column) |
| `summa_openWQ_allocspace.f90` | allocates the start-of-timestep prognostic snapshot (`gru_hru_dom_doubleVec`) OpenWQ needs for flux mass balance |
| `openWQ.f90`, `openWQInterface.f90` | Fortran <-> C++ (`iso_c_binding`) wrapper around the OpenWQ class |
| `OpenWQ_hydrolink.{cpp,h}`, `OpenWQ_interface.{cpp,h}` | C++ hydrolink that drives OpenWQ's coupler calls |
| `CMakeLists.txt` | builds the `openWQ` object library (OpenWQ C++ + hydrolink) and wires its dependencies |

`openwq/` (the OpenWQ C++ source) and `_deps/` (locally-built third-party libraries) are **git-ignored build inputs** — they are recreated with the steps below, not committed. `build/cmake_build_openwq/` (the CMake build tree) is git-ignored for the same reason as `build/cmake_build/`.

---

## Dependencies

OpenWQ's `develop` branch requires, in addition to SUMMA's normal deps (NetCDF, LAPACK/BLAS):

| dependency | notes |
|------------|-------|
| **OpenWQ**    | `ue-hydro/openwq`, branch `develop`, cloned into `build/source/openwq/openwq` |
| **Armadillo** | must be built with `ARMA_USE_HDF5` — OpenWQ headers use `hid_t` / `arma::hdf5_name`, pulled in transitively through `<armadillo>` only when that macro is set |
| **PhreeqcRM** | `usgs-coupled/phreeqcrm` — geochemistry engine, a hard dependency of `src/global/OpenWQ_wqconfig.hpp` |
| **SUNDIALS/CVODE** | OpenWQ's chemistry / sediment ODE solver (`src/compute/solver_sundials`) includes `<cvode/cvode.h>` unconditionally. SUMMA's own SUNDIALS build (IDA/KINSOL) already ships CVODE + nvecserial + core |
| **HDF5**      | OpenWQ output + external-forcing input. HDF5 >= 1.12 (incl. the 2.x series) needs `-DH5_USE_110_API`, which `CMakeLists.txt` adds to the `openWQ` target |

> OpenWQ upstream only officially supports Linux builds (see `openwq/containers/` for the
> Docker / Apptainer recipe). The steps below are a native macOS build (MacPorts toolchain:
> `/opt/local/bin/{gfortran,gcc,g++}`, MacPorts HDF5, Accelerate for LAPACK) and have been
> verified to compile, link, and run `summa_sundials_openwq.exe -v`.

---

## 1. Fetch OpenWQ and build the local dependencies

```sh
cd build/source/openwq
./build_deps.sh          # clones openwq/, builds Armadillo + PhreeqcRM into _deps/
```

`build_deps.sh` is idempotent. Override the compilers / install location with env vars
(`CXX`, `CC`, `DEPS_DIR`) — see the top of the script. It prints the exact `cmake` configure
line to use when it finishes.

To do it by hand instead:

```sh
cd build/source/openwq

# OpenWQ C++ source
git clone -b develop https://github.com/ue-hydro/openwq.git openwq

# Armadillo (with HDF5 support)
mkdir -p _deps && cd _deps
curl -L -o armadillo.tar.xz \
  "https://sourceforge.net/projects/arma/files/armadillo-14.2.3.tar.xz/download"
tar -xf armadillo.tar.xz && cd armadillo-14.2.3
# enable ARMA_USE_HDF5 in the config template so the wrapper library is built with it
perl -0pi -e 's{^// #define ARMA_USE_HDF5$}{#define ARMA_USE_HDF5}m' \
  include/armadillo_bits/config.hpp.cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/local/bin/g++ -DCMAKE_C_COMPILER=/opt/local/bin/gcc \
  -DCMAKE_PREFIX_PATH=/opt/local \
  -DCMAKE_CXX_FLAGS="-DH5_USE_110_API -I/opt/local/include" \
  -DCMAKE_SHARED_LINKER_FLAGS="-L/opt/local/lib -lhdf5" \
  -DCMAKE_INSTALL_PREFIX="$PWD/../armadillo-install"
cmake --build build -j4 && cmake --install build
cd ..

# PhreeqcRM
git clone --depth 1 https://github.com/usgs-coupled/phreeqcrm.git phreeqcrm && cd phreeqcrm
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/local/bin/g++ -DCMAKE_C_COMPILER=/opt/local/bin/gcc \
  -DBUILD_SHARED_LIBS=ON -DPHREEQCRM_BUILD_TESTS=OFF -DPHREEQCRM_FORTRAN_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/../phreeqcrm-install"
cmake --build build -j4 && cmake --install build
```

---

## 2. Configure and build SUMMA with OpenWQ

`-DUSE_OPENWQ=ON` produces `summa_openwq.exe`, or `summa_sundials_openwq.exe` when combined
with `-DUSE_SUNDIALS=ON` (works with both). From `build/`:

```sh
DEPS=$PWD/source/openwq/_deps
SUNDIALS=<path to your SUNDIALS install>     # e.g. .../SummaSundials/sundials/instdir

FC=gfortran LIBRARY_LINKS="-framework Accelerate" \
cmake -B cmake_build_openwq -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_SUNDIALS=ON \
  -DUSE_OPENWQ=ON \
  -DSPECIFY_LAPACK_LINKS=ON \
  -DSUNDIALS_DIR="$SUNDIALS/lib/cmake/sundials" \
  -DPhreeqcRM_DIR="$DEPS/phreeqcrm-install/lib/cmake/PhreeqcRM" \
  -DCMAKE_PREFIX_PATH="/opt/local;$DEPS/armadillo-install;$DEPS/phreeqcrm-install;$SUNDIALS"

cmake --build cmake_build_openwq -j4
```

The executable is written to `bin/summa_sundials_openwq.exe`. CMake bakes the `_deps/*/lib`
directories into the executable's RPATH, so it runs without setting `DYLD_LIBRARY_PATH`:

```sh
./bin/summa_sundials_openwq.exe -v
```

To actually run the coupled model, OpenWQ additionally needs its JSON configuration and a
file named `openwq_mainJSONFile_fullPath.txt` (containing the full path to the OpenWQ master
JSON) in the working directory where the executable is launched.

---

## Keeping `CMakeLists.txt` in sync with OpenWQ

`CMakeLists.txt` in this directory mirrors the source-file list and the dependency wiring
(PhreeqcRM, SUNDIALS/CVODE, HDF5) from OpenWQ's own `openwq/CMakeLists.txt` (the
`summa_openwq` target). When OpenWQ is updated, re-check that list against
`openwq/CMakeLists.txt`.
