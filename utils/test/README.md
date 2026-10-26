# SUMMA tests

Tests and example domains for SUMMA, grouped by what they exercise:

- [`test_regression/`](test_regression/README.md) -- correctness and performance
  regression: does a code change alter SUMMA's answers (stable vs. development,
  serial vs. MPI, mizuRoute on vs. off), and how does MPI scale. Requires
  private reference data and multiple local builds; not runnable from this
  repository alone.
- [`test_mizuroute/`](test_mizuroute/README.md) -- everything about the
  SUMMA-mizuRoute coupling: a synthetic-network smoke test, a real-network
  test compared against t-route, and a real-world calibration/evaluation
  case (Bow River at Banff).
- [`test_calibration/`](test_calibration/README.md) -- parameter calibration
  (`summa*_opt.exe`, needs `-DUSE_MPI=ON -DUSE_MIZUROUTE=ON`):
  `test_calibration_bow.sh`, a short DDS calibration of the Bow at Banff in
  2004 (inside the GRACE record, so streamflow and storage targets can be
  scored together); `test_calibration_bow_pareto.sh`, the same case scored on
  discharge and GRACE by NSGA-II and by DDS on one budget; `test_calibration_gulkana_wolverine.sh` and
  `test_calibration_wolverine.sh`, NSGA-II on two glacier basins, one GRU and reach each,
  scoring discharge, GRACE, mass balance and a synthetic well; and `multi_case_example/`, which calibrates several
  cases from one manifest on a stub dataset built from the bundled Bow domain.
- [`test_ngen/`](test_ngen/readme.md) -- example NextGen case studies (Provo,
  gauge_01073000) showing how a SUMMA setup looks under NextGen, run either
  standalone or coupled through ngen with t-route routing. Not a pass/fail
  test suite; `test_mizuroute/` borrows its bundled domains as test inputs.

Runnable right after cloning, given a build with the options named above
(see [docs/index.md](../../docs/index.md)): `test_mizuroute/test_mizuroute_bundled.sh`,
`test_calibration/test_calibration_bow.sh`, `test_calibration/test_calibration_bow_pareto.sh`,
`test_calibration/multi_case_example/make_stub_century.bash`
followed by the run it prints.
