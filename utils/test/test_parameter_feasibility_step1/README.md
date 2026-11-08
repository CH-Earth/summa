# Generic calibration feasibility: production regression tests

These tests build the production modules in `build/source/calibration/` and exercise
nested `hspatial(:)%dat(:)` physical fields, ordered constraints, spatial multipliers,
parameter ratios, and SUMMA-specific field extraction and write-back.

From the repository root:

```sh
./utils/test/test_parameter_feasibility_step1/run_tests.sh
```

Override compiler and flags with `FC=gfortran` and `FFLAGS='-O2 -g -fcheck=all'`.
The runner uses a temporary build directory, so the checkout is not polluted by
`.o` or `.mod` files. The SUMMA adapter test uses actual SUMMA data and mapping
types with a stub only for unrelated `globalData` state; it is **not** an MPI or
full SUMMA integration test.

## Production contract

- Call `initialize_summa_calibration_adapter` after SUMMA initializes its native
  parameters. The generic constraint plan checks fixed-only ordered chains once.
- Call `is_feasible(context,x,feasible,err,message)` for each candidate, where
  `x` contains untransformed calibration variables. An ordinary rejection returns
  `feasible=.false., err=0`; a setup/checker failure returns `err/=0`.
- Random/DDS and NSGA-II retain their own existing retry loops. They call the
  same predicate without knowing what physical constraint failed.
- The worker applies **accepted** candidates using `apply_summa_calibration_adapter`
  without repeating calibration ordered checks. Structural and finite-value
  safeguards remain on the model application path.
- This refactor intentionally removes scalar bound propagation. Search bounds
  stay at configured values; search trajectories are not bit-for-bit invariant.

`parameter_fields` constructs nested physical fields, `parameter_constraints`
evaluates their validity, `model_calibration` exposes the generic checker, and
`summa_parameter_fields` only handles SUMMA storage/mapping.

The historical Step 1 stand-alone source copies and block-based parity script
were removed after production integration; the previous prototype remains
available in Git history (commit `d6dc3ca6`).
