# What's new

This page provides simple, high-level documentation about what has changed in each new release of SUMMA.

## Version 3.0.4 (pre-release)
- Initial addition of the "What's new" page
- Added pull request template
- Adds HRU/GRU info to error messages
- Fixes a segfault of mysterious origin when using JRDN snow layering
- Fixes a water balance error w.r.t transpiration
- Fixes the output message to report the correct solution type
- Adds tolerance to balance check in updatState.f90
- Changes all float data types to `rk`, for "real kind", which is intended to make it easier to switch from double to single precision
- Computes wall clock time for each HRU and time step
- Fixes a bug that incorrectly writes scalarTotalET and scalarNetRadiation to output in cases where canopy calculations are skipped
- Canopy ice content check in check_icond.f90 now generates a warning if ice > 0 for T > 0 instead of a graceful exit. Graceful exit still exists if ice > 1E-3.
- Added case_study folder and Reynolds Mountain East albedo decay experiment
- Fixes a bug where solar angle incorrectly gets set to 0 during polar days
- Add deflate (compression level) option to outputControl file -- default level is 4 if not specified
- Fixes an unnecessary rounding error on SAI and LAI values in PHENOLOGY routine
- Fixes a bug where the SUMMA version is incorrectly reported by "summa.exe -v"
- Fixes a bug that incorrectly writes scalarRainPlusMelt to output in cases where snow layers do not exist
- Changed part "(a,1x,i0)" to "(a,1x,i0,a,f5.3,a,f5.3)" in check_icond.f90 line 277 to print out error correctly.
- Fixes a bug where restart files are not read correctly in cases where the parallelization argument `-g` for setups that have >1 HRU per GRU
