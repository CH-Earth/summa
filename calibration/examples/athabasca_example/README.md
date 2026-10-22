# Example: calibration and strong scaling of the Athabasca basin (69 HRUs, 15 parameters)

Scripts used to calibrate the semi-distributed Athabasca domain (Athabasca_sedi_3: 69 HRUs = 69 river
segments, outlet segment 222) with asynchronous parallel DDS, calling SUMMA + coupled mizuRoute (KW)
through `summa_python_interface.J(params)`, and to measure strong scaling on Purdue Anvil.

The calibration engine is a custom asynchronous DDS (Tolson & Shoemaker 2007) in an mpi4py
master-worker driver: rank 0 proposes candidates, ranks 1..N-1 each run one `J()` evaluation at a
time in a child process (`eval_one*.py`), so a fatal model error ends only that evaluation.

For comparison, the same `J()` setup is also driven by SPOTPY (SCE-UA, DDS, LHS) in MPI mode
(`scripts/spotpy/`); only the optimizer differs. Results of both are below.

## Files

| File | Purpose |
|---|---|
| `scripts/calibrate_athabasca_asyncdds_mpi.py` | Async DDS driver, summa-spotpy build (the wrapper in this PR) |
| `scripts/calibrate_athabasca_asyncdds_mpi_develop.py` | Same driver for the NCAR/summa `develop` port, time-series output off |
| `scripts/calibrate_athabasca_asyncdds_mpi_scaling_develop_v2.py` | Scaling driver: all workers start at once; only successful SUMMA runs count toward the budget; canopy-invalid candidates are redrawn |
| `scripts/eval_one.py`, `scripts/eval_one_develop.py` | One `J()` evaluation in its own process (old build / develop) |
| `scripts/spotpy/calibrate_athabasca_spotpy.py` | SPOTPY driver (`sceua`, `dds`, `lhs`, `parallel="mpi"`) on the same `J()`, worker setup and failure handling as the async DDS develop driver |
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
