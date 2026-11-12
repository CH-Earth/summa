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

## Oct 30-Nov 3: why does the Fortran DDS stop at 0.806? (v2 with its search rules)

- The Fortran DDS differs from v2 in two rules: it starts from a random set (v2: defaults) and rejects
  out-of-bounds candidates (v2: reflects). Each was switched on in v2 alone and together
  (seed 42, 5,000 runs, 49 workers).

| Run | Job | Start | Bounds | @1000 | @2000 | @3000 | Best |
|---|---|---|---|---|---|---|---|
| v2 | 21608093 | defaults | reflect | 0.841 | 0.875 | 0.884 | 0.8925 |
| def_reject | 21895156 | defaults | reject | 0.781 | 0.785 | 0.819 | 0.8698 |
| rand_reject (both rules) | 21895099 | random | reject | 0.816 | 0.854 | 0.855 | 0.8564 |
| rand_reflect | 21895100 | random | reflect | 0.740 | 0.744 | 0.744 | 0.7813 |
| Fortran DDS | 21621404 | random | reject | 0.802 | 0.805 | 0.805 | 0.8056 |

- With both rules, v2 still reached 0.856, so the two rules alone do not explain 0.806. The Fortran best
  sits in a different parameter region (opposite ends of qSurfScale and routingGammaScale, low
  frozenPrecipMultip) from the v2 runs.

## Nov 2-9: Pymoo GA, four crossover/mutation settings

- Driver options added to `scripts/pymoo/calibrate_athabasca_pymoo.py`: `--cx-prob`, `--cx-eta`, `--mut-prob`,
  `--mut-eta` (SBX crossover and polynomial mutation). Same J(), 5,000 runs, seed 42, start = defaults + Latin
  hypercube sets, one whole node per job, workers = offspring per generation.

| | S1 baseline | S2 more diversity | S3 more exploration | S4 balanced | GA Oct |
|---|---|---|---|---|---|
| Job | 21896879 | 21896880 | 21896881 | 21896882 | 21595996 |
| Population / generations | 40 / 125 | 100 / 50 | 50 / 100 | 50 / 100 | 98 / 101 |
| SBX prob / eta | 0.9 / 15 | 0.9 / 15 | 0.9 / 5 | 0.9 / 10 | 0.9 / 15 |
| PM prob / eta | 0.1 / 20 | 0.1 / 20 | 0.2 / 10 | 0.15 / 15 | 1/15 / 20 |
| Best KGE | 0.837 | 0.853 | 0.880 | 0.797 | 0.882 |
| Idle workers | 15% | 29% | 9% | 6% | 22% |

- S3 and S4 start from the same population, so the 0.880 vs 0.797 difference comes from the operators alone.
  S1 and S4 lost most of their population diversity within 1,000-1,500 runs and stalled (S4 in a separate
  region: qSurfScale ~48, routingGammaScale ~14); S3 kept the most diversity and kept improving.
- No setting beat the October GA (0.882) or async DDS v2 (0.8925). The October GA (large population) and S3
  (strong operators) reach ~0.88 by different routes.
- All idle time comes from waiting at generation boundaries; replaying the same run times asynchronously gives
  ~1% idle at 43 evaluations per worker, as for DDS.
- v2, the October GA and S3 reach 0.88-0.89 with very different parameter sets (k_soil 1.2e-6 to 9.5e-3,
  heightCanopyTop 6.5-100 m); routingGammaScale is at its upper bound in all three. About half of the runs
  with KGE < -1 have frozenPrecipMultip > 1.9, in every scenario.

## Nov 3-8: Fortran DDS with other seeds

- Seed patch in `summa_parameter_sampling.f90` (`scripts/fortran_dds/patches/fortran_dds_seed.patch`): the seed is
  read from `SUMMA_DDS_SEED` (default 42) and printed in the log; nothing else changed. The rebuilt `libsumma.so` is byte-identical to the original in
  code and data (only the build time stamp and Anvil's link metadata differ).
- Seeds 43-45 (jobs 21987864-66; the first submission failed at start-up because `coldState.nc` was not
  copied into the new state folders). The jobs stopped at their 30 h limit before 5,000 runs:

| Seed | Runs | @1000 | @2000 | @3000 | Best |
|---|---|---|---|---|---|
| 42 (original) | 5,077 | 0.802 | 0.805 | 0.805 | 0.806 |
| 43 | 3,983 | 0.760 | 0.821 | 0.828 | 0.849 |
| 44 | 4,933 | 0.745 | 0.757 | 0.769 | 0.823 |
| 45 | 4,370 | 0.822 | 0.827 | 0.829 | 0.847 |

- With other seeds the same code reaches 0.82-0.85, so 0.806 was partly a bad draw. These runs went at
  133-164 runs/h, against 239 runs/h for the original; the cause is not known yet.
- The program writes `_calibration.nc` only at a normal end, so runs stopped by the time limit keep their
  KGE values (in the log) but not their parameter sets. Writing the file during the run (`nf90_sync`) would
  avoid this.

## Nov 9-10: Fortran DDS started from the SUMMA defaults

- Option added to the patched program (`scripts/fortran_dds/patches/fortran_dds_seed_defstart.patch`, separate
  executable `summa_mizuroute_opt_defstart.exe`, job script `slurm/run_fortran_dds_defstart.slurm`):
  `SUMMA_DDS_START=default` starts the search from the default values instead of a random set and prints
  the 15 starting values. Search rules and everything else unchanged.
- Jobs 22039052-54 (seeds 42-44), all finished normally: 5,000 runs, no failures.

| Seed | Start | @1000 | @2000 | @3000 | @5000 | Runs/h | Median run |
|---|---|---|---|---|---|---|---|
| 42 | random (original) | 0.802 | 0.805 | 0.805 | 0.806 | 239 | 11.6 min |
| 42 | defaults | 0.775 | 0.780 | 0.804 | 0.867 | 214 | 12.9 min |
| 43 | defaults | 0.825 | 0.856 | 0.883 | 0.902 | 222 | 12.2 min |
| 44 | defaults | 0.792 | 0.810 | 0.814 | 0.891 | 216 | 12.1 min |

- From the defaults, the unmodified search reaches 0.867-0.902 (v2: 0.8925). With seed 42, the only change
  is the start: 0.806 -> 0.867. The random start is the main reason for the low score.
- Most of the gain comes after 3,000 runs, and all three found their best within the last 100 runs, so they
  had not converged at 5,000. The default start is not faster per run.
- Same-seed reruns of the random-start jobs (43b-45b, jobs 22038122-24) are still running. They measure the
  spread from asynchronous worker timing (so far 43b is at 0.767 where seed 43 reached 0.849 at the same
  stage) and will give complete 5,000-run random-start results.

## Nov 10-12: final results, best parameter sets and run speed

- The same-seed reruns 43b-45b (jobs 22038122-24) finished normally: 5,000 runs each, no failures. Their
  first model run has the same KGE as the first attempts, so they share seed and starting point; they
  replace the time-limited seeds 43-45 below. Script: `scripts/fortran_dds/analyze_fortran_dds.py`.

| Seed | Start | Runs | @1000 | @3000 | Best | Wall time | Median run |
|---|---|---|---|---|---|---|---|
| 42 (original) | random | 5,077 | 0.802 | 0.805 | 0.806 | 20.8 h | 11.6 min |
| 43 (43b) | random | 5,100 | 0.695 | 0.767 | 0.769 | 43.7 h | 25.6 min |
| 44 (44b) | random | 5,133 | 0.754 | 0.824 | 0.864 | 42.0 h | 24.8 min |
| 45 (45b) | random | 5,094 | 0.828 | 0.831 | 0.834 | 36.3 h | 21.2 min |
| 42 | defaults | 5,113 | 0.775 | 0.804 | 0.867 | 23.3 h | 12.9 min |
| 43 | defaults | 5,135 | 0.825 | 0.883 | 0.902 | 22.5 h | 12.2 min |
| 44 | defaults | 5,140 | 0.792 | 0.814 | 0.892 | 23.2 h | 12.1 min |

- Every default-start run ended above every random-start run (worst 0.867 vs best 0.864). With 4
  independent random starts against 3 default starts this is suggestive (about 1 in 35 by chance), not
  conclusive. The random start remains the most likely cause of the 0.806.
- Seed 43 from the defaults (0.902) is the best Athabasca result so far. It reached KGE 0.88 in 2,886 runs
  and 0.85 in 8.5 h (v2: 2,973 runs and 10.3 h).
- Same seed, same start and same code still ended up to 0.08 apart (seed 43: 0.849 in the time-limited
  run, 0.769 in the rerun): asynchronous worker timing changes the search path. Setups need several seeds
  each to be compared.
- Best parameter sets: the three default-start runs agree (qSurfScale 2.7-3.7, routingGammaScale ~99,000,
  routingGammaShape 3.9-4.4, frozenPrecipMultip 0.94-1.36, aquiferBaseflowRate ~1e-7, heightCanopyTop
  2-4 m). The random-start runs end in two other places: qSurfScale high with routingGammaScale low
  (seed 42 original, 43b), or frozenPrecipMultip 0.61-0.75 with a baseflow rate 10-100 times higher (44b,
  45b). routingGammaScale is at its upper bound (1e5) in 5 of the 7 best sets.
- Run speed is set by heightCanopyTop, not by the job or the executable (all builds are Release -O3):
  runs with heightCanopyTop >= 20 m take 18-26 min against 11-13 min below, in every job. The 15
  parameters explain 97.5% of the variance of log run time (random forest, cross-validated), 95% of it
  from heightCanopyTop. The reruns spent 99% of their runs at 50-97 m, which explains the slow seeds 43-45.
- Scheduling is not the problem: in all seven jobs, 0.9-1.3% of the time comes after the last dispatch
  and 0.4-0.8% of worker capacity is idle.

## Next

- Update this branch to current develop and bring in the develop version of `summa_c_api.f90`.
- Repeat the develop calibration with 3-5 seeds (async DDS and SPOTPY DDS); raise the theta_sat lower
  bound to ~0.39 and widen the routingGammaScale upper bound.
- If SCE-UA is repeated: fewer complexes than workers, or a larger budget.
- Repeat v2 and v1 with 3-5 seeds, and the GA with S3-type settings (2-3 seeds).
- Cap frozenPrecipMultip at about 2 (about half of all KGE < -1 runs lie above 1.9).
- Rerun the seven best sets with time-series output for a KGE breakdown (r, alpha, beta).
- Cap heightCanopyTop at about 30 m (runs above 20 m take twice as long; the best sets use 2-4 m).
- Suggest to the SUMMA develop Fortran DDS: an option to start from the defaults, and writing
  `_calibration.nc` during the run. (Develop's Nov 9 update already sets the seed through `run_index`.)
- Validate the best sets outside 2016-2019, starting with seed 43 from the defaults (0.902).
