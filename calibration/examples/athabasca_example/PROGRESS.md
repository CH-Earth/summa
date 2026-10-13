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

## Next

- Update this branch to current develop and bring in the develop version of `summa_c_api.f90`.
- Repeat the develop calibration with 3-5 seeds; raise the theta_sat lower bound to ~0.39.
