# What's new
This page provides simple, high-level documentation about what has changed in each new release of SUMMA. Please add any changes made in pull requests to under the `Pre-release` header. Use `Minor changes` sub-heading for changes that do not affect science outputs or are likely to affect only a minority of users. Use `Major changes` for anything else.

A commit-level list of every change on the 4.x line is kept in
[`docs/assets/changes_fromV3Summa.txt`](assets/changes_fromV3Summa.txt); the summary below only
covers the user-facing highlights.

## Version 4.0.0 (experimental)

### Build system
- SUMMA is now built with CMake. Options select the SUNDIALS solvers (`-DUSE_SUNDIALS=ON`),
  the NextGen framework (`-DUSE_NEXTGEN=ON`), the OpenWQ water-quality coupling
  (`-DUSE_OPENWQ=ON`), and the build type (`-DCMAKE_BUILD_TYPE=Release|Debug`). See the
  [installation instructions](installation/SUMMA_installation.md).

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
- The number of soil layers no longer has to be the same in every HRU.

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
- Optional coupling to the OpenWQ water-quality framework (`build/source/openwq/`).
- Runs as a NextGen submodule; NextGen test cases are in `test_ngen/`.
- Large refactor: object-oriented flux routines, much shorter `computFlux.f90` and the
  individual flux modules, simplified Jacobian assembly.

## Pre-release
### Major changes
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
