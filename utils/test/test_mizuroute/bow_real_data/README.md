# SUMMA-mizuRoute coupled test case: Bow River at Banff (CAN_05BB001)

A real, lumped-catchment SUMMA-mizuRoute coupling case, distinct from the
[toy problem and Provo real-network tests](../README.md) one level up: it
uses a real (single-HRU lumped) river network and remapping, ~9 years of real
distributed forcing, and observed daily streamflow for objective-function
evaluation -- rather than either a synthetic network or a network converted
purely to compare against t-route.

## Layout
- `summa_inputs/` -- forcing, attributes, parameters, decisions for the Bow
  at Banff domain. Two forcing sets, each with its own list:
  `CAN_05BB001_em_earth_distributed_1980-1990.nc` (`forcingFileList.txt`),
  hourly EM-Earth on Mountain Standard Time, declared as `-7:00` in its time
  units so that `tmZoneInfo = ncTime` places the sun correctly, and
  `CAN_05BB001_rdrs_2002-2017.nc` (`forcingFileList_rdrs.txt`), hourly RDRS
  v2.1 on UTC, 2002-04 to 2017-12, which overlaps GRACE. The RDRS file joins
  the `Bow_at_Banff_multivar` basin-average forcing, with `hruId` set to 
  this domain's 101 and `data_step` added.
- `observations/` -- observations for the RDRS period, both on UTC: daily
  flow (`CAN_05BB001_daily_flow_2002-2017.nc`) and the JPL GRACE anomaly
  (`CAN_05BB001_grace_tws.nc`). The basin sits inside
  one 3-degree mascon, so the cell holding its centroid is its basin average.
- `mizuroute_inputs/` -- the real river network (`topology.nc`), the
  lumped-catchment-to-routing-HRU remapping (`lumped_to_hru.nc`), mizuRoute
  parameters (`mizuroute.param`), and observed daily flow for the EM-Earth
  period (`CAN_05BB001_daily_flow_observations.nc`, on Mountain Standard Time
  to match that forcing).
- `settings/summa_fileManager.txt` -- SUMMA file manager for this domain.
- `settings/mizu_control_CAN_05BB001.toml` -- mizuRoute coupling
  configuration: hydrofabric, remapping, observations, and the `[objective]`
  section (KGE against observed flow over 1982-10-01 to 1983-10-01).

## Running

Requires a build configured with `-DUSE_MIZUROUTE=ON` (see
[docs/index.md](../../../../docs/index.md)), run from the repository root
(the file manager and TOML config use paths relative to it):

```
bin/summa_sundials_mizuroute.exe -m utils/test/test_mizuroute/bow_real_data/settings/summa_fileManager.txt \
                                  -c utils/test/test_mizuroute/bow_real_data/settings/mizu_control_CAN_05BB001.toml
```

Output goes to `work/` (created if needed, gitignored). With `write_aligned`
set in the TOML, mizuRoute also writes an aligned evaluation time series and
the objective-function value alongside the routed output.
