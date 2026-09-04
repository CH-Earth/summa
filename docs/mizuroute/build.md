# Building mizuRoute

mizuRoute can be used with SUMMA in either **coupled** or **standalone**
mode. The build procedure depends on the intended workflow.

## Git submodules

mizuRoute and other external dependencies are included in the SUMMA repository
as Git submodules. After cloning SUMMA, initialize the submodules with:

```bash
git submodule update --init --recursive
```

This retrieves the versions of mizuRoute and the other submodules associated
with the current SUMMA revision.


## Coupled SUMMA--mizuRoute

For coupled simulations, mizuRoute is compiled directly as part of the SUMMA
CMake build. No separate mizuRoute build step is required.

Enable the coupling with:

```bash
-DUSE_MIZUROUTE=ON
```

For example, on macOS using Homebrew:

```bash
export FC=/opt/homebrew/bin/gfortran
export LIBRARY_LINKS='-llapack'
export SUNDIALS_DIR=../../../sundials/build/

cmake -B ../cmake_build -S ../. \
    -DUSE_MPI=ON \
    -DUSE_SUNDIALS=ON \
    -DUSE_MIZUROUTE=ON \
    -DSPECIFY_LAPACK_LINKS=ON \
    -DCMAKE_BUILD_TYPE=Debug

cmake --build ../cmake_build --target all -j
```

When `USE_MIZUROUTE` is enabled, the SUMMA build system compiles the required
mizuRoute source and the SUMMA--mizuRoute coupling layer and links them into
the SUMMA library. The mizuRoute source is included in the SUMMA repository
as a Git submodule.

The other CMake options shown above are application dependent; only
`USE_MIZUROUTE=ON` is required specifically to enable mizuRoute coupling.

mizuRoute support is optional. When `USE_MIZUROUTE` is not enabled, SUMMA is
built without the coupled routing components and can be run using the standard
SUMMA workflow without a mizuRoute dependency.

## Standalone mizuRoute

For a sequential workflow, mizuRoute is built as a standalone executable.
The authoritative build instructions are provided in the
[mizuRoute documentation](https://mizuroute.readthedocs.io/en/main/users_guide/Build_model.html).

SUMMA also includes an example build script:

```text
utils/build/make_mizuRoute.sh
```

The script provides an example of building the mizuRoute submodule and its
ParallelIO dependency on macOS using Homebrew. It builds the native mizuRoute
`route_runoff` executable and copies the resulting executable to the SUMMA
`bin` directory.

Because the standalone build relies on the native mizuRoute build system,
users should refer to the mizuRoute documentation for platform-specific
compiler and dependency requirements.
