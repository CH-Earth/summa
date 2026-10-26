# gulkana_wolverine — two glacier basins, two GRUs, two reaches

Gulkana and Wolverine glaciers, Alaska, as one SUMMA domain of 2 GRUs, each draining to its own
mizuRoute reach. It is the domain of `../test_calibration_gulkana_wolverine.sh` and
`../test_calibration_wolverine.sh`, which test that a calibration runs, scores and ranks trials on
several targets, GRUs and reaches; they calibrate nothing.

| | |
| --- | --- |
| Forcing | ERA5, Oct 2015 – Jul 2018, 34 monthly files |
| GRU 1, Gulkana | HRUs 1–3 (26.6 km²), draining to reach 1 |
| GRU 2, Wolverine | HRUs 5–7 (24.7 km²), draining to reach 2 |

| Observation file | Record | Scored by |
| --- | --- | --- |
| `gulkana_USGS_15478040_daily_flow.nc` | USGS 15478040 discharge, m3 s-1 | both GRUs' test, reach 1 |
| `wolverine_USGS_15236900_daily_flow.nc` | USGS 15236900 discharge, open-water season | both tests, reach 2 |
| `gulkana_grace_tws.nc`, `wolverine_grace_tws.nc` | GRACE JPL mascon anomaly at each GRU's centroid, mm | both GRUs' test |
| `wolverine_mass_balance_seasonal.nc` | glacier-wide seasonal balance, mm w.e. | Wolverine test, GRU 2 |
| `wolverine_well_hru5_daily.nc` | synthetic well: HRU 5's aquifer storage in a run at the default parameters, m | Wolverine test, HRU 5 |

The calibration spins up over the year before `start_time`, from the cold state, so the tests score
October 2016 to July 2018. A GRACE anomaly is referenced to the scored period (`baseline_start`,
`baseline_end`), since the run cannot reproduce the product's own baseline.

The reaches are estimates.
