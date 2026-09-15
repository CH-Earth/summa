# SUMMA tests

Tests and example domains for SUMMA, grouped by what they exercise:

- [`test_regression/`](test_regression/README.md) -- correctness and performance
  regression: does a code change alter SUMMA's answers (stable vs. development,
  serial vs. MPI, mizuRoute on vs. off), and how does MPI scale. Requires
  private reference data and multiple local builds; not runnable from this
  repository alone.
- [`test_mizuroute/`](test_mizuroute/README.md) -- SUMMA-mizuRoute coupling
  tests. `test_mizuroute_bundled.sh` is fully self-contained (bundled domain,
  synthetic river network) and can run from a fresh checkout. The other
  scripts route the real Provo hydrofabric network, to compare against the
  same domain routed through ngen/t-route.
- [`test_ngen/`](test_ngen/readme.md) -- example NextGen case studies (Provo,
  gauge_01073000) showing how a SUMMA setup looks under NextGen, run either
  standalone or coupled through ngen with t-route routing. Not a pass/fail
  test suite; `test_mizuroute/` borrows its bundled domains as test inputs.

Only `test_mizuroute/test_mizuroute_bundled.sh` is expected to run unmodified
right after cloning the repo (it needs a build configured with
`-DUSE_MIZUROUTE=ON`, see [docs/index.md](../../docs/index.md)).
