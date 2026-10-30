# Example: calibration and strong scaling of the Athabasca basin (69 HRUs, 15 parameters)

Scripts used to calibrate the semi-distributed Athabasca domain (Athabasca_sedi_3: 69 HRUs = 69 river
segments, outlet segment 222) with asynchronous parallel DDS, calling SUMMA + coupled mizuRoute (KW)
through `summa_python_interface.J(params)`, and to measure strong scaling on Purdue Anvil.

The calibration engine is a custom asynchronous DDS (Tolson & Shoemaker 2007) in an mpi4py
master-worker driver: rank 0 proposes candidates, ranks 1..N-1 each run one `J()` evaluation at a
time in a child process (`eval_one*.py`), so a fatal model error ends only that evaluation.

For comparison, the same `J()` setup is also driven by SPOTPY (SCE-UA, DDS, LHS) in MPI mode
(`scripts/spotpy/`) and by a Pymoo genetic algorithm (`scripts/pymoo/`); only the optimizer differs.
Async DDS v2 adds the search rules of SUMMA develop's own Fortran DDS, and that Fortran DDS
(`summa_mizuroute_opt.exe`, used unmodified) was also run on the same basin. Results of all are below.

## Files

| File | Purpose |
|---|---|
| `scripts/calibrate_athabasca_asyncdds_mpi.py` | Async DDS driver, summa-spotpy build (the wrapper in this PR) |
| `scripts/calibrate_athabasca_asyncdds_mpi_develop.py` | Same driver for the NCAR/summa `develop` port, time-series output off |
| `scripts/calibrate_athabasca_asyncdds_mpi_scaling_develop_v2.py` | Scaling driver: all workers start at once; only successful SUMMA runs count toward the budget; canopy-invalid candidates are redrawn |
| `scripts/eval_one.py`, `scripts/eval_one_develop.py` | One `J()` evaluation in its own process (old build / develop) |
| `scripts/spotpy/calibrate_athabasca_spotpy.py` | SPOTPY driver (`sceua`, `dds`, `lhs`, `parallel="mpi"`) on the same `J()`, worker setup and failure handling as the async DDS develop driver |
| `scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py` | Async DDS v2: v1 plus log-space search, ordered constraints checked before a run (redrawn, no run spent), the default set in a full first wave; options `--bounds reflect\|reject`, `--start defaults\|random` |
| `scripts/athabasca_dds_v2_settings.toml` | v2 search settings in the develop calibration TOML format: log k_soil, canopy and soil ordered constraints |
| `scripts/pymoo/calibrate_athabasca_pymoo.py` | Pymoo GA (SBX + polynomial mutation, tournament selection; population 98, 49 offspring per generation) on the same `J()`, with a checkpoint after every generation for `--resume` |
| `scripts/fortran_dds/setup_fortran_dds.py` | Writes a Fortran DDS run folder (`fortran_dds/<tag>/`: `fileManager.txt`, `config.toml`) and `settings/SUMMA/localParamInfo_fortran.txt` from the J() config and the v2 settings, so both search the same space |
| `scripts/fortran_dds/convert_fortran_dds_output.py` | Converts the Fortran DDS `*calibration.nc` into the `evaluations.csv` / `summary.json` files the Python drivers write |
| `scripts/load_balance_athabasca.py` | Drain fraction and idle-worker fraction of finished runs |
| `scripts/plot_worker_timeline_athabasca.py` | One panel per run: each model run as a bar on its worker, coloured by KGE |
| `scripts/plot_load_balance_scaling.py` | Measured run times in blocks of 1,000 runs, replayed with asynchronous scheduling on 192/95/47/23 workers |
| `scripts/compare_athabasca_optimizers.py` | Compares finished runs (async DDS, SPOTPY, ...): table of best KGE, wall time, runs per hour, runs/hours to reach given KGE levels, plus a best-so-far figure vs runs and vs hours |
| `scripts/athabasca_params.csv` | The 15 calibration parameters, defaults and bounds |
| `scripts/make_athabasca_obs_nc.py` | Converts the observed daily streamflow CSV to NetCDF |
| `scripts/run_best_hydrograph_develop.py` | Full-period run (2014-09 to 2024-12) of one parameter set with output on |
| `scripts/postprocess_hydro_develop.py` | Reads `Q_reach` at segment 222 from SUMMA's `*_timestep.nc`, daily means, KGE for 2016-2019 and 2020-2024, hydrograph |
| `scripts/plot_scaling.py` | Strong-scaling figure (wall time and parallel efficiency vs nodes) |
| `scripts/develop_port/validate_basin.py` | Same parameters on old build vs develop (output on / off) |
| `scripts/develop_port/test_handshake.py` | Python -> C -> Fortran handshake test (`summa_echo`) |
| `settings/*.example.toml` | J() configs: summa-spotpy format, and develop format with `write_timeseries` true / false |
| `slurm/run_athabasca_spotpy.slurm` | SPOTPY runs: `<algo> <repetitions> [ngs] [smoke]` (usage examples in the file header) |
| `slurm/run_athabasca_asyncdds_develop_v2.slurm`, `slurm/run_athabasca_pymoo.slurm` | Async DDS v2 and Pymoo GA runs (usage in the file headers) |
| `slurm/build_develop_opt.sh` | Builds the Fortran DDS (`summa_mizuroute_opt.exe`) and a single-run executable from the same develop checkout as `libsumma.so`, in a separate build folder |
| `slurm/run_fortran_dds.slurm`, `slurm/run_fortran_dds_check.slurm`, `slurm/run_crosscheck.slurm` | Fortran DDS run; one run with the defaults (must give J()'s KGE 0.5528); best sets cross-checked under both spin-up methods |
| `slurm/*.slurm`, `slurm/submit_scaling.sh`, `slurm/build_develop.sh` | Anvil job scripts and the develop build |

## Before running

Replace the placeholders with your own values:

- `<WORK_ROOT>`: folder that holds `summa-spotpy/` (this repository, with `my_work_athabasca/`) and `summa-develop/`
- `<SCRATCH>`: scratch space for per-worker SUMMA output (scaling runs)
- `<ACCOUNT>`: Slurm allocation

The scripts expect the basin folder `my_work_athabasca/` with `settings/`, `observations/`, `scripts/`
(these files), `results/` and `output/`. Forcing, observations and results are not included in the repository.

`summa_python_interface.py` now accepts an optional `SUMMA_LIB` environment variable to load another
`libsumma` (e.g. the develop port); without it the default library is used as before.

## Notes

- Observations: use daily means stamped at the **end** of each 24 h period (period-end file), which is
  how `J()` labels simulated days. A start-stamped file shifts the comparison by one day.
- Coupled mizuRoute in `develop` supports KW, Muskingum-Cunge and diffusive wave only; IRF stops with
  "not implemented: use standalone mizuRoute".
- In `develop`, routed flow is written as `Q_reach(time, seg, method)` inside SUMMA's `*_timestep.nc`.
- About 4-5% of evaluations fail: `theta_sat` below ~0.383 (SUMMA paramCheck) and canopy top <= bottom.

## Results (Anvil, Oct 2026)

| Run | Evaluations | Wall time | Best KGE 2016-2019 | KGE 2020-2024 |
|---|---|---|---|---|
| summa-spotpy build, 49 workers | 5,000 | 36.0 h | 0.8927 (0.8904 on period-end obs) | 0.903 |
| develop, output off, 49 workers | 5,000 | 28.4 h | 0.8806 | 0.895 |
| SYMFLUENCE best parameters run through `J()` | - | - | 0.820 | 0.804 |

Strong scaling (develop, output off, 127 workers per node, 5,000 successful evaluations):

| Nodes | Wall time | Speed-up | Parallel efficiency |
|---|---|---|---|
| 1 | 12.33 h | 1.00x | 100% |
| 2 | 6.34 h | 1.94x | 97% |
| 4 | 3.67 h | 3.36x | 84% |
| 8 | 2.26 h | 5.46x | 68% |

## Optimizer comparison on the same J() (develop, output off, 49 workers on 1 node, 5,000-run budget, seed 42)

| | Async DDS (job 21217337) | SPOTPY 1.6.7 DDS (job 21343300) | SPOTPY 1.6.7 SCE-UA, 49 complexes (job 21339837) |
|---|---|---|---|
| **Best KGE 2016-2019** | **0.8806** | 0.8510 | 0.7239 |
| Wall time | **28.4 h** | 32.4 h | 95.7 h |
| Model runs executed | 5,000 | 5,000 | 8,636 (5,000 kept) |
| Runs per hour | **176** | 154 | 90 |
| Failed runs | 4.2% | 3.6% | 8.7% |
| Best KGE after 1,000 / 2,500 runs | 0.790 / 0.857 | 0.801 / 0.836 | 0.647 / 0.698 |
| Best KGE after 12 h / 24 h | 0.841 / 0.880 | 0.816 / 0.846 | 0.698 / 0.698 |
| Hours to KGE 0.80 / 0.85 / 0.8748 | 10.5 / 14.2 / 22.7 | 7.4 / 30.1 / - | - / - / - |

Made with `scripts/compare_athabasca_optimizers.py` (every executed run counts, in the order it finished).

- Both DDS runs use the same rule (Tolson & Shoemaker 2007, r = 0.2, start from the defaults) and are
  asynchronous; async DDS ran 14% more model runs per hour and ended higher. SPOTPY DDS led until about
  1,600 runs while async DDS sat at 0.79. With one seed each, the 0.03 KGE gap is not yet evidence that
  one DDS implementation searches better; repeat runs with 3-5 seeds would show the spread.
- SCE-UA with 49 complexes (one per worker) spends 1,519 runs on burn-in and fits only about two shuffling
  loops into 5,000 runs; SPOTPY also completes each loop past the budget (3,636 runs executed and
  discarded). This setup is a poor fit for a 5,000-run budget rather than a verdict on SCE-UA.
- The best DDS sets agree on the runoff regime: low qSurfScale (2.6-2.9) and routingGammaScale at its upper
  bound (about 1e5) in both, so that bound should be widened in the next runs.

## Six calibrations compared (Oct 22-28; develop, output off, 49 workers on 1 node, 5,000-run budget, seed 42)

| | Async DDS v1 | Async DDS v2 | Fortran DDS | Pymoo GA | SPOTPY DDS | SPOTPY SCE-UA |
|---|---|---|---|---|---|---|
| Job | 21217337 | 21608093 | 21621404 | 21595996 | 21343300 | 21339837 |
| **Best KGE 2016-2019** | 0.8806 | **0.8925** | 0.8056 | 0.8821 | 0.8510 | 0.7239 |
| Model runs | 5,000 | 5,000 | 5,000 | 4,998 | 5,000 | 8,636 |
| Wall time | 28.4 h | 41.5 h | **20.9 h** | 35.6 h | 32.4 h | 95.7 h |
| Runs per hour | 176 | 120 | **239** | 140 | 154 | 90 |
| Failed runs | 4.2% | 3.4% | 0.0% | 1.8% | 3.6% | 8.7% |
| Runs to KGE 0.80 / 0.85 / 0.88 | 1,621 / 2,355 / 4,131 | 555 / **1,255** / **2,973** | 893 / - / - | **442** / 1,970 / 4,335 | 927 / 4,630 / - | - / - / - |
| Hours to KGE 0.85 | 14.2 | **10.3** | - | 18.9 | 30.1 | - |
| Best KGE after 1,000 / 2,500 runs | 0.790 / 0.857 | 0.841 / 0.876 | 0.802 / 0.805 | 0.820 / 0.859 | 0.801 / 0.836 | 0.647 / 0.698 |

Made with `scripts/compare_athabasca_optimizers.py`; the Fortran run was first converted with
`scripts/fortran_dds/convert_fortran_dds_output.py`. Every executed run counts, in the order it finished.

- **Async DDS v2** ends highest and needs the fewest runs to reach 0.85 and 0.88. Compared with v1 it takes
  log-space k_soil, the canopy check before a run and the full first wave from the Fortran DDS, and keeps v1's
  start from the defaults and reflection at the bounds. It is the slowest DDS per hour because its runs took longer.
- **Fortran DDS** (SUMMA develop's `summa_mizuroute_opt.exe`, unmodified): the default set gives KGE
  0.55281681983953113, the same as `J()` to 16 digits. It ran the most runs per hour with no failures, partly
  because all trials share one spin-up year, but stopped improving at about 0.805 after ~1,000 runs (random
  first wave, out-of-bound steps rejected).
- **Pymoo GA** reaches 0.80 first (442 runs) and ends at 0.8821; its workers wait at each generation boundary.
- One seed each: the top three are within 0.012 KGE, so the ranking is suggestive, not proven.
- Fortran DDS setup notes: the Fortran driver stops the whole job on the first failed SUMMA run, so
  theta_sat's lower bound is raised to 0.39 (below ~0.383 SUMMA fails); the state folder and spin-up start are
  set per run folder; `seg_outlet` is removed, because in subset mode every model run writes a topology file
  into the shared `hfabric_path` (without it the largest-area reach, 222, is scored).
- Several best values sit at a bound (routingGammaScale at its upper bound in four runs; tempCritRain and
  Fcapil near their lower bounds), and similar KGE comes from quite different parameter sets (v1 vs v2
  k_soil 1.9e-4 vs 1.2e-6 m/s). The bounds should be revisited before repeat runs.
