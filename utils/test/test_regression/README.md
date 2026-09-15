# SUMMA code regression tests

These scripts check that a SUMMA code change has not altered model output --
"regression" in the sense of comparing runs, not a bundled pass/fail test
suite. All of them read paths under `$HOME/data/...` and `$HOME/models/...`
that are **not part of this repository** (private reference domains and
side-by-side builds); edit the "User settings" block at the top of each
script before running it.

## What each script does

- `parallel_sims.sh` / `parallel_sims_original.sh` -- run the same domain with
  a stable build, a development build (serial), and a development build under
  MPI at several process counts (`parallel_sims_original.sh` is the same
  comparison against the pre-refactor MPI implementation). Writes one output
  file per run/decomposition plus a log per run.
- `regression_test.sh` -- compares the output files those runs produced:
  mizuRoute on vs. off, stable serial vs. development serial, development
  serial vs. MPI np=1, and development serial vs. each decomposed MPI run
  (reassembling the per-process HRU/GRU ranges to diff against the matching
  slice of the serial run). Reports the maximum absolute difference per HRU-
  and GRU-level variable and an overall PASS/FAIL per comparison; PASS means
  bit-identical.
- `test_mizuroute_coupling.sh` -- compares SUMMA's coupled mizuRoute output
  (`q_reach`) against a standalone mizuRoute run (`KWroutedRunoff`) on the
  same real river network, reach by reach. For a version of this comparison
  that runs entirely from the repository (synthetic network, no external
  data or standalone mizuRoute build), see
  [`../test_mizuroute/test_mizuroute_bundled.sh`](../test_mizuroute/test_mizuroute_bundled.sh).
- `parallel_time.sh` -- summarizes wall-clock time, speedup and parallel
  efficiency across the MPI process counts `parallel_sims.sh` ran, comparing
  the current MPI implementation against a saved "original" run.
- `plotvars.R` -- plots a handful of state variables at one HRU from two of
  the output files above, for visually inspecting where two runs diverge.

## Requirements

`ncdiff`, `ncap2`, `ncks`, `ncwa`, `ncrename` (NCO), `bc`, and for
`parallel_sims.sh`, `mpirun` and MPI-enabled SUMMA builds. `plotvars.R` needs
the R `ncdf4` package.
