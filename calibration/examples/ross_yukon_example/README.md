# Example: calibration of the Ross basin, Yukon (15 parameters)

Scripts used for a 1000-evaluation DDS calibration of the Ross basin with SPOTPY,
calling SUMMA + mizuRoute through `summa_python_interface.J(params)`.

| File | Purpose |
|---|---|
| `scripts/calibrate_ross_spotpy.py` | DDS calibration driver (SPOTPY) |
| `scripts/ross_params.csv` | The 15 calibration parameters and their bounds |
| `scripts/test_J_ross.py`, `scripts/validate_two_points.py` | Validation of `J()` on Ross |
| `scripts/export_iteration_results.py`, `scripts/plot_iterations.py`, `scripts/check_progress.py` | Progress monitoring and results |
| `scripts/make_ross_obs_nc.py` | Converts the observed streamflow CSV to NetCDF |
| `settings/ross_config.example.toml` | Config template; replace `<ROSS_DIR>` with your local folder |

The scripts locate files relative to this folder (`settings/`, `observations/`, `results/`).
The Ross forcing, observations and results are not included in the repository.
