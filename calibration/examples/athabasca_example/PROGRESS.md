# Athabasca example: progress log

Work since the last push to PR #651 (Oct 1, 2026). Newest last. Each step matches one commit on
`fhan/objfunc-c-api-python`.

## Oct 6-8: first full Athabasca calibration (summa-spotpy build)

- Async DDS driver with mpi4py (rank 0 = DDS, ranks 1..N-1 = SUMMA), each `J()` call in its own
  process so a fatal model error only fails that evaluation.
- 15 parameters, KGE on 2016-2019 at outlet segment 222, coupled mizuRoute KW, 49 workers.
- Job 21150407: 5,000 evaluations in 36.0 h (139 evaluations/h); best KGE 0.8927.
  SYMFLUENCE (same basin) stopped at 3,200 evaluations after 47.9 h with 0.8748.
- About 5% of runs fail: theta_sat below ~0.383 (paramCheck) and canopy top <= bottom.

## Oct 8-9: observation time stamps

- Observed daily means are now stamped at the end of each 24 h (period-end file), as `J()` labels
  simulated days. Same flows, dates one day later. On this file the job 21150407 best set scores 0.8904.

## Oct 9: port of the C wrapper to NCAR/summa develop

- Reason: `write_timeseries = false` (no time-series output during calibration) exists only in develop.
- `summa_c_api.f90` ported to develop's new `evaluate_objective` signature (on Anvil; to be brought
  into this branch together with the update to current develop), plus a `summa_echo` test function.
- `summa_python_interface.py`: optional `SUMMA_LIB` to load another `libsumma`; default unchanged.
- `test_handshake.py`: all 15 values reach Fortran in order at full precision; checksums match.

## Oct 9: develop port validated (same parameters, old build vs develop)

| Basin | old | develop, output on | develop, output off |
|---|---|---|---|
| Athabasca | 0.8904 (20.5 min) | 0.8904 (14.8 min) | 0.8904 (13.9 min, 0 MB output) |
| Ross | 0.9133 | 0.9127 | 0.9127 |

- develop renamed scalarSnowSublimation to scalarGroundSublimation; `validate_basin.py` writes a
  converted output list.
- IRF routing is not available in coupled develop mizuRoute (KW, MC and DW only), so KW is kept.

## Oct 10-11: full calibration on develop, output off

- Job 21217337, same settings as 21150407: 5,000 evaluations in 28.4 h (176 evaluations/h, 21% less
  time); best KGE 0.8806 (period-end observations).
- Different search path (plateau at 0.79 until 1,650 evaluations), same final region
  (low qSurfScale, slow hillslope routing at the routingGammaScale upper bound).

## Oct 12: full-period runs 2014-2024 (job 21304052)

- `run_best_hydrograph_develop.py` + `postprocess_hydro_develop.py` (develop writes `Q_reach` into
  SUMMA's `*_timestep.nc`).

| Parameter set | KGE 2016-2019 | KGE 2020-2024 |
|---|---|---|
| summa-spotpy build best | 0.890 | 0.903 |
| develop best | 0.881 | 0.895 |
| SYMFLUENCE best, run through J() | 0.820 | 0.804 |

## Oct 10-13: strong scaling on whole Anvil nodes (develop, output off)

- Pilot on one node: per-run time rises 1.8x from 16 to 127 workers, but a full node gives the most
  throughput. Quick failures used up to 51% of the budget at start-up with 127 workers.
- v2 scaling driver: canopy-invalid candidates redrawn, failed SUMMA runs recorded but not counted,
  so every node count does exactly 5,000 successful runs (check run: 256 counted, 9 failed, 4 redrawn).

| Nodes | Workers | Wall time | Speed-up | Efficiency |
|---|---|---|---|---|
| 1 | 127 | 12.33 h | 1.00x | 100% |
| 2 | 255 | 6.34 h | 1.94x | 97% |
| 4 | 511 | 3.67 h | 3.36x | 84% |
| 8 | 1,023 | 2.26 h | 5.46x | 68% |

- Loss at 8 nodes: early-search evaluations take ~25-27 min vs ~16.5 min later, and with 1,023
  workers more of the budget is spent before the search improves (~22 points); idle workers at the
  end of the run (~13 points).

## Oct 13-14: SPOTPY driver on the same J()

- `scripts/spotpy/calibrate_athabasca_spotpy.py`: SPOTPY 1.6.7 `sceua`, `dds` and `lhs` with
  `parallel="mpi"`, reusing the develop driver's per-rank setup, child-process `J()` call, 2 h timeout,
  -1e6 penalty and canopy check. DDS: r = 0.2, `x_initial` = defaults.
- Seeding: SPOTPY reseeds numpy when the sampler is created, so the seed is passed as
  `random_state = 42 + rank`. SPOTPY's DDS evaluates the starting set on rank 0, so rank 0 gets the
  same SUMMA setup as the workers (found in a local MPI test with a stand-in model).
- Smoke tests on Anvil passed for SCE-UA and DDS. The full DDS run's first evaluation (the default
  parameters) gives KGE 0.5528, the same as in job 21217337, so both drivers run the same J() setup.

## Oct 13-17: full SPOTPY runs (1 shared node, 49 workers, 5,000-run budget, seed 42)

- DDS, job 21343300: 5,000 runs in 32.4 h (154 runs/h); best KGE 0.8510.
- SCE-UA with 49 complexes, job 21339837: 8,636 runs executed (SPOTPY keeps the first 5,000) in 95.7 h
  (90 runs/h); best KGE 0.7239 (0.7202 among the 5,000 kept).

## Oct 22: comparison with async DDS (`scripts/compare_athabasca_optimizers.py`)

| | Async DDS | SPOTPY DDS | SPOTPY SCE-UA |
|---|---|---|---|
| Best KGE 2016-2019 | 0.8806 | 0.8510 | 0.7239 |
| Wall time | 28.4 h | 32.4 h | 95.7 h |
| Runs per hour | 176 | 154 | 90 |
| Runs / hours to KGE 0.85 | 2,355 / 14.2 h | 4,630 / 30.1 h | not reached |

- Async DDS: fastest, highest final KGE, and the only run past SYMFLUENCE's 0.8748 (after 22.7 h).
- SPOTPY DDS reached 0.80 first (927 runs vs 1,621) and then fell behind after about 1,600 runs; one seed
  each, so search-path luck and implementation cannot be separated yet.
- SCE-UA: burn-in of 1,519 runs and about two shuffling loops within the budget; not a fair test of
  SCE-UA at this budget.

## Oct 22-24: Pymoo genetic algorithm on the same J()

- `scripts/pymoo/calibrate_athabasca_pymoo.py`: rank 0 runs the GA and hands out sets to free workers;
  population 98 (defaults + 97 Latin hypercube sets), 49 offspring per generation, 4,998 runs in total.
- Smoke test 21591237 passed; full run 21595996 (whole node): best KGE 0.8821 in 35.6 h. None of the 97
  random starting sets had a positive KGE; 79 of the 91 failures were canopy-invalid sets with no SUMMA run.

## Oct 23-26: async DDS v2 (job 21608093)

- `scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py` with `scripts/athabasca_dds_v2_settings.toml`:
  log-space k_soil (1e-8 to 1e-2), ordered constraints checked before a run, defaults + 48 DDS steps in the
  first wave. Bounds reflected and start from the defaults, as in v1.
- 5,000 runs in 41.5 h: best KGE 0.8925 (v1 0.8806); KGE 0.85 after 1,255 runs (v1 2,355); idle workers
  0.4% (v1 2.2%).

## Oct 23-26: SUMMA develop's own Fortran DDS on Athabasca

- Code used unmodified: `summa_driver_opt.f90`, `summa_parameter_sampling.f90`, `parameter_search.f90`,
  built as `summa_mizuroute_opt.exe` from the same develop checkout as `libsumma.so` (`slurm/build_develop_opt.sh`).
- Inputs from the J() setup (`scripts/fortran_dds/setup_fortran_dds.py`): same 15 parameters, bounds, 5,000
  runs, log k_soil, canopy constraint.
- Check 21604889: the default set gives KGE 0.55281681983953113, identical to `J()` to 16 digits.
- Smoke tests found four setup fixes: spin-up start, state folder, theta_sat lower bound 0.39 (a failed run
  stops the whole job), and no `seg_outlet` (in subset mode each run writes a topology file into the shared
  `hfabric_path`).
- Full run 21621404: 5,000 runs in 20.9 h (239 runs/h), no failures, best KGE 0.8056. It stops improving
  after about 1,000 runs.

## Oct 26-28: six-way comparison and load balancing

- Table and best-so-far figure from `scripts/compare_athabasca_optimizers.py` (README, "Six calibrations
  compared"). v2 is best (0.8925), then GA 0.8821 and v1 0.8806.
- `scripts/plot_worker_timeline_athabasca.py`: DDS runs keep every worker busy; the GA shows idle stripes at
  generation boundaries; SCE-UA scores KGE <= -1 in most runs and leaves most workers idle twice for hours.
- `scripts/load_balance_athabasca.py`, `scripts/plot_load_balance_scaling.py`: replaying 1,000-run blocks
  on 192/95/47/23 workers gives an idle fraction of about 14% at 5 evaluations per worker, falling to about
  1% at 43.

## Next

- Update this branch to current develop and bring in the develop version of `summa_c_api.f90`.
- Repeat the develop calibration with 3-5 seeds (async DDS and SPOTPY DDS); raise the theta_sat lower
  bound to ~0.39 and widen the routingGammaScale upper bound.
- If SCE-UA is repeated: fewer complexes than workers, or a larger budget.
- Repeat v2, GA and v1 with 3-5 seeds; run the Fortran DDS from the defaults with reflection at the bounds
  to separate the code from these two choices.
- Validate the best sets outside 2016-2019.
