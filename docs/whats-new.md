# What's new

This page provides simple, high-level documentation about what has changed in each new release of SUMMA.

Documented here are large changes that will likely affect the majority of users:
- New features
- Changes to the functionality of the code
- Changes to default settings
- Changes in code conventions (e.g. variable renaming)
- Etc.

## Version 3.0.4 (pre-release)

- Initial addition of the "What's new" page
- Added pull request template
- Adds HRU/GRU info to error messages
- Fixes a segfault of mysterious origin when using JRDN snow layering
- Fixes a water balance error w.r.t transpiration
- Fixes the output message to report the correct solution type
- Adds tolerance to balance check in updatState.f90
- Changes all float data types to `rk`, for "real kind", which is intended to
  make it easier to switch from double to single precision.
