# What's new
This page provides simple, high-level documentation about what has changed in each new release of SUMMA. Please add any changes made in pull requests to under the `Pre-release` header. Use `Minor changes` sub-heading for changes that do not affect science outputs or are likely to affect only a minority of users. Use `Major changes` for anything else.

A commit-level list of every change on the 4.x line is kept in
[`docs/assets/changes_fromV3Summa.txt`](assets/changes_fromV3Summa.txt); the summary below only
covers the user-facing highlights.

## Version 4.0.0 (experimental)

### Build system
- SUMMA is now built with CMake. Options select the SUNDIALS solvers (`-DUSE_SUNDIALS=ON`),
  the NextGen framework (`-DUSE_NEXTGEN=ON`), the OpenWQ water-quality coupling
  (`-DUSE_OPENWQ=ON`), MPI (`-DUSE_MPI=ON`), mizuRoute (`-DUSE_MIZUROUTE=ON`), and the build
  type (`-DCMAKE_BUILD_TYPE=Release|Debug`). See the
  [installation instructions](installation/SUMMA_installation.md).
- Executable names record the options they were built with (`summa_sundials_mizuroute.exe`,
  and so on), so differently configured builds can share `bin/` without overwriting one
  another. `CMAKE_BUILD_TYPE` is not part of the name, so Debug and Release builds of the
  same options still collide; give them separate build directories. Each build script also
  takes a `clean` argument.

### Numerical solution
- New `num_method` options `kinsol` and `ida` use the SUNDIALS KINSOL and IDA solvers
  (require a SUNDIALS build). The built-in backward-Euler solver is `homegrown` (the legacy
  name `itertive` still works).
- `fDerivMeth` now selects a numerical vs. analytical Jacobian; the in-flux-routine numerical
  derivatives were removed. Analytical derivatives are required for `kinsol`/`ida`.
- New `nrgConserv` decision (renamed from `howHeatCap`): temperature/closed-form heat capacity
  vs. an enthalpy formulation with a soil temperature–enthalpy lookup table or analytical
  relation.
- New solver-control parameters for the backward-Euler step count and for SUNDIALS
  relative/absolute tolerances and IDA controls; all have sensible defaults.
- The moisture-based form of Richards' equation was removed — `f_Richards` only does the
  mixed form now.

### Spatial representation
- HRUs can be subdivided into multiple spatial **domains** (upland plus glacier
  accumulation / clean-ablation / debris-ablation columns; wetland and lake columns are
  scaffolded but not yet active). Data structures gained a `dom` dimension throughout; the
  restart and attributes files gained `dom`, `glac` and glacier-grid variables. Runs without
  glaciers or wetlands are unchanged.
- Glacier area evolves during a run: once a year a shallow-ice flow model with debris
  transport is run on a per-glacier grid and the glacier domains are rebuilt from the result
  (see [glacier geometry updates](configuration/SUMMA_configuration.md#glacier-geometry-updates)).
  Four new GRU parameters control the debris model; all have defaults.
- The number of soil layers no longer has to be the same in every HRU.

### Parallelization
- New horizontal (spatial) HRU/GRU domain decomposition: a build configured with
  `-DUSE_MPI=ON` distributes GRUs across MPI ranks, each writing its own range of GRUs/HRUs
  to output, alongside the existing serial executable (PR #631). Reassembled MPI output is
  bit-identical to serial output.
- Spatial indexing distinguishes three reference domains throughout: the *file* domain
  (`nGRU_file`), the *run* domain selected with `-g` (`startGRU_domain`, `nGRU_domain`), and
  the *local* domain assigned to a rank (`startGRU_local`, `nGRU_local`).
- A second level of parallelism runs concurrent calibration trials across ranks, alongside
  the domain decomposition within each trial (PR #639).

### Calibration and objective functions
- SUMMA can evaluate a streamflow objective function at the end of a run and write it to the
  output file: `kge`, `kgep`, `nse`, `mae` or `rmse`, with optional `log`, `power` or
  `box-cox` transformation of the flows, over a chosen evaluation period (PR #637). It is
  computed only when observations are configured, so ordinary runs are unaffected.
- Model parameters can be overridden from the command line with `--param <name> <value>`,
  without editing the trial-parameter file (PR #637).
- The simulation lifecycle is now reusable in-process (`initialize_summa` / `run_summa` /
  `finalize_summa`), so one executable can run many parameter sets without restarting
  (PR #637).
- New automatic calibration driver (`summa[_sundials]_opt.exe`) using DDS sampling, with
  parameter transformations, ordered-parameter constraints, a shared spinup restart, and
  NetCDF output of every trial (PR #639).
- Many basins can be calibrated in one job through a manifest file (`--manifest`) that lists
  cases and a configuration template (PR #639).

### Process options
- New `infRateMax` decision for the maximum infiltration rate (`topmodel_GA`, `GreenAmpt`,
  `noInfExc`), and new `surfRun_SE` decision for saturation-excess surface runoff
  (`homegrown_SE`, `FUSEPRMS`, `FUSEAVIC`, `FUSETOPM`, `zero_SE`). Distinct
  `scalarSurfaceRunoff_IE` / `scalarSurfaceRunoff_SE` output fluxes were added.
- New `aquiferIni` decision (`fullStart` / `emptyStart`).
- Wind-profile / stability changes to fix over-estimated snow sublimation (affects
  `veg_traits = CM_QJRMS1988` most).
- Soil and snow longwave emissivity updated (0.98/0.99 → 0.96/0.98).

### Input / output
- Simulation start/end time (`simStartTime`, `simEndTime`) and `tmZoneInfo` are set in the
  file manager, not the model decisions file. The file manager version string is
  `SUMMA_FILE_MANAGER_V3.0.0`.
- Configuration can now be given as a TOML file passed with `-c`, covering the paths and
  settings the file manager held plus the mizuRoute, observation and calibration settings.
  The legacy control file is still accepted (`-m`, now also spelled `--control`)
  (PRs #632, #639).
- New `read_force` decision (buffered vs. per-step forcing reads; was `readForcing`) and
  `write_buff` decision (buffered vs. per-step output writes; was `writeOutput`).
- Fluxes and soil compression are written as means over the output window rather than the
  value at the end of the last sub-step.
- New energy/mass balance output variables (`balanceCasNrg`, `balanceVegNrg`, `balanceSnowNrg`,
  `balanceSoilNrg`, `balanceVegMass`, `balanceSnowMass`, `balanceSoilMass`, `balanceAqMass`),
  `meanStepSize`, and the GRU-level `basin__StorageChange`.
- Routing-histogram variables are no longer written by default (set `allowRoutingOutput` to
  re-enable). Several derived heat-capacity/conductivity scalars that no longer point to
  anything were removed from the output list.

### Other
- Optional coupling to mizuRoute for river network routing of basin runoff
  (`-DUSE_MIZUROUTE=ON`), adding `Q_reach`, `q_basin` and `upArea` output. Coupled mizuRoute
  needs the whole domain on one process, so it rejects `-g` and MPI domain parallelism
  (PR #632).
- Optional coupling to the OpenWQ water-quality framework (`build/source/openwq/`).
- Runs as a NextGen submodule; NextGen test cases are in `utils/test/test_ngen/`.
- Large refactor: object-oriented flux routines, much shorter `computFlux.f90` and the
  individual flux modules, simplified Jacobian assembly.

## Pre-release
### Major changes
- Fixed: writing a restart file failed for any run with more than one GRU and a glacier grid.
- A soil layer's lateral outflow is capped at the drainable water it holds over the data step.
- The infiltration closure under saturation scales to the depth of the zone it averages over, so
  a thin glacier debris column closes off as a soil column does.
- Exfiltration returns the surplus a nearly full column cannot store whatever face it arrived
  through, not only the lateral inflow. This closes the water balance of a glacier debris column
  fed by melt from below; a free-draining soil column is unchanged, since the relief only opens
  once available storage falls below 2 mm.
- General cleanup and shortening of computFlux.f90, vegNrgFlux.f90, snowSoilNrgFlux.f90, vegLiqFlux.f90, snowLiqFlux.f90, soilLiqFlux.f90, groundwatr.f90, and bigAquifer.f90 
- Added object-oriented methods to simplify flux routine calls in computFlux and improve modularity
    - classes for each flux routine were added to data_types.f90
    - large associate statemements are no longer needed in computFlux (associate blocks are now much shorter)
    - the length of computFlux has been decreased substantially
- Added a new decision to set maximum infiltration rate method
- Bug fix: fixed a problem with snow sublimation due to a bug in transitioning from exponential to log wind profile below canopy.

### Minor changes
- Updated SWE balance check in coupled_em for cases where all snow melts in one of the substeps

## Version 3.2.0
### Major changes
- Addition to compute wall clock time for each HRU and time step
- Fixes a bug that incorrectly writes scalarTotalET and scalarNetRadiation to output in cases where canopy calculations are skipped
- Added case_study folder and Reynolds Mountain East albedo decay experiment
- Fixes a bug where restart files are not read correctly in cases with the parallelization argument `-g` for setups that have >1 HRU per GRU

### Minor changes
- Fixes a bug where solar angle incorrectly gets set to 0 during polar days
- Canopy ice content check in check_icond.f90 now generates a warning if ice > 0 for T > 0 instead of a graceful exit. Graceful exit still exists if ice > 1E-3.
- Add deflate (compression level) option to outputControl file -- default level is 4 if not specified
- Fixes an unnecessary rounding error on SAI and LAI values in PHENOLOGY routine
- Fixes a bug where the SUMMA version is incorrectly reported by "summa.exe -v"
- Fixes a bug that incorrectly writes scalarRainPlusMelt to output in cases where snow layers do not exist
- Changed part "(a,1x,i0)" to "(a,1x,i0,a,f5.3,a,f5.3)" in check_icond.f90 line 277 to print out error correctly.
- Adds scalarSnowDrainage variable when melting of the snow without a layer
- Changes the logic for creating the first snow layer: instead of creating the layer when snow-without-a-layer exceeds the maximum depth of the 1st layer, the first layer is now created if snow-without-a-layer exceeds the average of specified 1st layer minimum and maximum depth (zminLayer1 and zmaxLayer1_lower in localParamInfo.txt)
- Added documentation of lookup table provenance

## Version 3.1.0
- Initial addition of the "What's new" page
- Added pull request template
- Adds HRU/GRU info to error messages
- Fixes a segfault of mysterious origin when using JRDN snow layering
- Fixes a water balance error w.r.t transpiration
- Fixes the output message to report the correct solution type
- Adds tolerance to balance check in updatState.f90
- Changes all float data types to `rk`, for "real kind", which is intended to make it easier to switch from double to single precision
