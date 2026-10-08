# OpenWQ Integration

This directory contains the SUMMA-side coupling to [OpenWQ](https://github.com/ue-hydro/openwq):

| file | purpose |
|------|---------|
| `summa_openWQ.f90`            | passes the water volumes and the water fluxes of SUMMA, and of the internally coupled mizuRoute when it is active, to OpenWQ; walks the `gru -> hru -> dom` structures (every (HRU, domain) pair is one OpenWQ cell column) |
| `openWQ.f90`, `openWQInterface.f90` | Fortran <-> C++ (`iso_c_binding`) wrapper around the OpenWQ class |
| `OpenWQ_hydrolink.{cpp,h}`, `OpenWQ_interface.{cpp,h}` | C++ hydrolink that declares the compartments and drives OpenWQ's coupler calls |
| `CMakeLists.txt` | builds the `openWQ` object library (OpenWQ C++ + hydrolink) and wires its dependencies |

The river side of the coupling lives with the mizuRoute coupling code:
`../mizuroute/wq_exchange.f90` builds the map from SUMMA GRUs to reaches, and
`../mizuroute/network_routing.f90` records the water budget of every reach over the SUMMA time step.

## What OpenWQ sees

One OpenWQ instance carries the solutes of both models.

| compartment | cells | content |
|-------------|-------|---------|
| `SCALARCANOPYWAT`       | (HRU, domain) columns x 1                | canopy storage |
| `ILAYERVOLFRACWAT_SNOW` | columns x max snow layers                | snow layers |
| `RUNOFF`                | columns x 1                              | liquid water reaching the surface in the step (rain plus melt) |
| `ILAYERVOLFRACWAT_SOIL` | columns x max soil layers                | soil layers |
| `SCALARAQUIFER`         | columns x 1                              | aquifer (a basin-wide aquifer uses the first column of the GRU) |
| `RUNOFF_TO_STREAM`      | columns x 1                              | runoff held while SUMMA routes it within the GRU |
| `ILAYERVOLFRACWAT_LAKE` | columns x max lake layers                | lake layers, only if a domain has lake layers |
| `RIVER_NETWORK_REACHES` | reaches x 1                              | mizuRoute reaches, only when mizuRoute is coupled (`simulation.use_mizuroute = true`) |

Cell ids are `<hruId>_z<layer>` (`<hruId>_d<domain>_z<layer>` when the HRU has several domains) and
`<segId>` for the reaches.

External water fluxes: `PRECIP`, and `GLACIER_ICE_MELT` when a domain has glacier ice.
Flux-concentration exports: `scalarRunoffVol_m3`, `averageRoutedRunoff`, `scalarTotalRunoff`, and
`Qlocal_out` (reach outflow) with mizuRoute.
Dependencies for the kinetic expressions: `SM`, `Tair_K`, `Tsoil_K`, `SWrad_Wm2`, `cellArea_m2`, `T`.

The header of `summa_openWQ.f90` lists the water fluxes that carry solute.

### Land to river mapping and its checks

Each GRU delivers its routed runoff (`averageRoutedRunoff`, the same value SUMMA hands to mizuRoute) to the
reaches in proportion to the inflow each reach receives from that GRU per unit runoff. The map is built once at
initialization by passing unit runoff through mizuRoute's own `remap_runoff` (if remapping is on) and `basin2reach`,
so it follows whatever hydrofabric, remapping file and `hw_drain_point` the run uses. Each reach then mixes its
start-of-step storage, upstream inflow, lateral inflow and any water-management injection, and passes the outflow
(and any abstraction) on, using the volumes mizuRoute accumulated over its sub-steps. The coupler verifies this
mapping on every run, independently of the case:

| when | check | message |
|------|-------|---------|
| init | runoff elements that reach no river reach | `WARNING: water quality: N runoff element(s) deliver to no river reach` |
| init | river-network area of each GRU vs its SUMMA area (tolerance 1%) | `WARNING: OpenWQ river coupling: the river-network area of N GRU(s) differs ...` |
| init | several routing methods active (water quality follows the first) | `OpenWQ river coupling: several routing methods are active ...` |
| end | runoff routed by SUMMA vs lateral inflow received by the reaches, outlet outflow, end storage, and the summed reach budget residual `start + upstream + lateral + injection - abstraction - outflow - end` | `OpenWQ river coupling, water check over the run (m3): ...` |

The lateral inflow should equal the routed runoff (0.00 %) and the residual should be at round-off (1e-15 of the
lateral inflow, verified on the Bow test with the kinematic wave, Muskingum-Cunge and diffusive wave methods
and remapping, and on the Athabasca case with 69 GRUs delivering to 69 reaches without remapping).
A residual far above that means mizuRoute moved water the coupler does not track (for example lake reaches with
precipitation and evaporation, which are not handled).

At initialization the coupler also writes `<RESULTS_FOLDERPATH>/openwq_compartments.json` (compartment names,
indices and cell dimensions, plus the river and lake compartment indices): the calibration scripts read it to
address land and river cells separately, since their cell indices overlap.

Tested paths: one HRU per GRU with routing (Bow), 69 GRUs with routing (Athabasca), and two HRUs in one GRU with
TOPMODEL lateral flow from the upslope HRU into the outlet HRU (`downHRUindex`): the solute follows the water into
the downslope soil, only the outlet HRU feeds the stream pool, and the tracer closes to 1e-8 kg (Bow `work/multihru`
and `work/multihru_noDown`).

Known limits: snow cells follow the layer index, not the snow layer itself, when SUMMA splits or merges layers;
solute moves one reach per time step; a reach shares the dependency values of the land column with the same
index when there are more columns than one; water management is handled (injection = solute-free water,
abstraction = solute leaving) but cannot be exercised through the internal coupling, whose TOML has no water
management options; mizuRoute lake reaches are not handled; lake, wetland and glacier domains have no test case
in the SUMMA tree and remain unverified.

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
