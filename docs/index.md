# Structure for Unifying Multiple Modeling Alternatives: SUMMA

[![Build Status](https://travis-ci.org/NCAR/summa.svg?branch=develop)](https://travis-ci.org/NCAR/summa)
[![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/NCAR/SUMMA/master/COPYING)
[![Documentation Status](https://readthedocs.org/projects/summa/badge/?version=latest)](http://summa.readthedocs.org/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.800772.svg)](https://doi.org/10.5281/zenodo.800772)

SUMMA (Clark et al., [2015a](#clark_2015a);[b](#clark_2015b);[c](#clark_2015c); [2021](#clark_2021)) is a hydrologic modeling framework that can be used for the systematic analysis of alternative model conceptualizations with respect to flux parameterizations, spatial configurations, and numerical solution techniques. It can be used to configure a wide range of hydrological model alternatives and we anticipate that systematic model analysis will help researchers and practitioners understand reasons for inter-model differences in model behavior. When applied across a large sample of catchments, SUMMA may provide insights in the dominance of different physical processes and regional variability in the suitability of different modeling approaches. An important application of SUMMA is selecting specific physics options to reproduce the behavior of existing models – these applications of "**model mimicry**" can be used to define reference (benchmark) cases in structured model comparison experiments, and can help diagnose weaknesses of individual models in different hydroclimatic regimes.

SUMMA is built on a common set of conservation equations and a common numerical solver, which together constitute the  “**structural core**” of the model. Different modeling approaches can then be implemented within the structural core, enabling a controlled and systematic analysis of alternative modeling options, and providing insight for future model development.

The important modeling features are:

 1. The formulation of the conservation model equations is cleanly separated from their numerical solution;

 2. Different model representations of physical processes (in particular, different flux parameterizations) can be used within a common set of conservation equations; and

 3. The physical processes can be organized in different spatial configurations, including model elements of different shape and connectivity (e.g., nested multi-scale grids and HRUs).


## Documentation
SUMMA documentation is available [online](http://summa.readthedocs.io/) and remains a work in progress. Additional SUMMA information including publications, test data sets, and sample applications can be found on the [SUMMA web site](http://www.ral.ucar.edu/projects/summa) at NCAR.


## Building SUMMA

SUMMA depends on NetCDF and LAPACK. Optional features pull in further dependencies, described below.

### Getting the source

Optional features live in git submodules under `external/`, so clone recursively:

```bash
git clone --recurse-submodules <repository-url>
```

If the repository is already cloned, fetch them with:

```bash
git submodule update --init --recursive
```

Each submodule is only needed by the option that uses it, and CMake fails at configure
time if you enable that option without it. A build with all options off needs no
submodules at all.

| Submodule | Needed by |
| --- | --- |
| [`parallel-utils`](https://github.com/CH-Earth/parallel-utils) | `USE_MPI=ON` |
| [`mizuRoute`](https://github.com/ESCOMP/mizuRoute) | `USE_MIZUROUTE=ON` |
| [`toml-f`](https://github.com/toml-f/toml-f) | `USE_MIZUROUTE=ON` (reads the TOML configuration file) |

To fetch just one:

```bash
git submodule update --init external/parallel-utils
```

### Configuring and building

Build with CMake from the `build` directory:

```bash
cd build
cmake -B cmake_build -S . -DUSE_SUNDIALS=ON
make -C cmake_build -j4
```

Executables are written to `bin/`.

Ready-made scripts for common platforms are in `build/cmake` (`build.mac.bash`, `build.cluster.bash`, `build.pc.bash`, and the `build_ngen.*` variants for NextGen). Run them from that directory, for example `./build.mac.bash`. Each sets the compiler and library paths for its platform; edit the options listed in the `cmake` call to change what gets built.

### Debug builds

All configurations default to `-DCMAKE_BUILD_TYPE=Release`. Change it to `Debug` for a build with `-Og`, backtraces, and array bounds checking:

```bash
cmake -B cmake_build -S . -DUSE_SUNDIALS=ON -DCMAKE_BUILD_TYPE=Debug
```

This applies to every configuration, NextGen included. Bounds checking makes the model considerably slower but turns out-of-range array access into an immediate, located error rather than silent memory corruption, so it is worth using when a run crashes or produces implausible values.

### Build options

Each is `OFF` by default and enabled with `-DOPTION=ON`.

| Option | Effect |
| --- | --- |
| `USE_SUNDIALS` | Build with the IDA and KINSOL solvers from the SUNDIALS suite. Required to use `num_method` of `ida` or `kinsol`. Needs `SUNDIALS_DIR` set to the SUNDIALS cmake directory if it is not on the default search path. |
| `USE_MPI` | Additionally build an MPI executable that distributes GRUs across ranks. Requires an MPI Fortran compiler and the `parallel-utils` submodule. Has no effect together with `USE_NEXTGEN`, which builds a library rather than an executable. |
| `USE_NEXTGEN` | Build the BMI library for the NextGen framework instead of the standalone executables. |
| `USE_OPENWQ` | Build with the OpenWQ water-quality coupler. |
| `USE_MIZUROUTE` | Build with mizuRoute river network routing. Requires the `mizuRoute` and `toml-f` submodules, and a TOML configuration file at run time. |
| `SPECIFY_LAPACK_LINKS` | Take LAPACK link flags from the `LIBRARY_LINKS` environment variable instead of detecting them automatically. |

The executable name records the options selected:

| Options | Executable |
| --- | --- |
| none | `summa.exe` |
| `USE_SUNDIALS` | `summa_sundials.exe` |
| `USE_OPENWQ` | `summa_openwq.exe` |
| `USE_SUNDIALS` + `USE_OPENWQ` | `summa_sundials_openwq.exe` |

`USE_MPI` does not replace the serial executable; it adds a second one alongside it, named `summa_mpi.exe` or `summa_sundials_mpi.exe`. MPI code is confined to a separate driver, so the serial executable carries no MPI dependency.


## Running SUMMA

Serial runs take the file manager as the only required argument:

```bash
./bin/summa_sundials.exe -m /path/to/fileManager.txt
```

Useful options are `-s <suffix>` to tag the output file names, `-g <startGRU> <countGRU>` to run a contiguous block of GRUs, and `-h <HRU>` to run a single HRU. Run the executable with no arguments for the full list.

### Running with MPI

The MPI executable takes the same arguments and partitions the run domain's GRUs evenly across ranks:

```bash
mpirun -np 4 ./bin/summa_sundials_mpi.exe -m /path/to/fileManager.txt
```

Each rank reads only its own GRUs and writes its own output file, tagged with the GRU range it covers, for example `run_1_G01-14_timestep.nc`. Reconstruct the full domain by concatenating the rank files along the `gru` and `hru` dimensions; all other dimensions are file-wide and identical across ranks, so no padding is needed.

Results are independent of the number of ranks. Serial output and reassembled MPI output agree bit-for-bit.

### Running with mizuRoute

A build configured with `USE_MIZUROUTE=ON` routes basin runoff through a river network and
writes reach-level streamflow into the usual SUMMA output files. It needs a TOML
configuration file in addition to the file manager, passed with `-c`:

```bash
./bin/summa_sundials.exe -m /path/to/fileManager.txt -c /path/to/config.toml
```

The TOML file carries the mizuRoute settings; the file manager continues to describe the
SUMMA side of the run. See `docs/mizuroute/` for the configuration format, the coupling
design, and the routing options.

The `-c` option only exists in builds configured with `USE_MIZUROUTE=ON`. Passing it to a
build without mizuRoute is an error rather than a silent no-op, so a run cannot quietly
ignore the configuration you gave it.

mizuRoute routing is currently serial: it is not combined with `USE_MPI`.


## Credits
SUMMA's initial implementation is described in two papers published in [Water Resources Research](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1944-7973). If you use SUMMA, please credit these two publications.

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, 2015a: A unified approach for process-based hydrologic modeling: Part 1. Modeling concept. _Water Resources Research_, [doi:10.1002/2015WR017198](http://dx.doi.org/10.1002/2015WR017198).<a id="clark_2015a"></a>

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015b: A unified approach for process-based hydrologic modeling: Part 2. Model implementation and case studies. _Water Resources Research_, [doi:10.1002/2015WR017200](http://dx.doi.org/10.1002/2015WR017200).<a id="clark_2015b"></a>
 
 * Clark, M. P., Zolfaghari, R., Green, K. R., Trim, S., Knoben, W. J. M., Bennett, A., Nijssen, B., Ireson, A., Spiteri, R. J., 2021: The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests. _Journal of Hydrometeorology_, [doi:10.1175/JHM-D-20-0175.1](http://dx.doi.org/10.1175/JHM-D-20-0175.1).<a id="clark_2021"></a>

In addition, an NCAR technical note describes the SUMMA implementation in detail:

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015c: The structure for unifying multiple modeling alternatives (SUMMA), Version 1.0: Technical Description. _NCAR Technical Note NCAR/TN-514+STR_, 50 pp., [doi:10.5065/D6WQ01TD](http://dx.doi.org/10.5065/D6WQ01TD).<a id="clark_2015c"></a>


## License
SUMMA is distributed under the GNU Public License Version 3. For details see the file `COPYING` in the SUMMA root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
