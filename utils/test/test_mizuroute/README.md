# SUMMA-mizuRoute coupling tests

Tests for the built-in SUMMA-mizuRoute coupling (`-DUSE_MIZUROUTE=ON`, see
[docs/index.md](../../../docs/index.md)). All three scripts run coupled SUMMA
end to end and check the routed output; they differ in which river network
they route through.

## `test_mizuroute_bundled.sh` -- synthetic network, fully self-contained

```
./test_mizuroute_bundled.sh [work_dir]
```

Runs the bundled `../test_ngen/domain_provo` domain with a synthetic,
single-chain river network built by `make_test_topology.py` (GRU *i* drains
to reach *i*, reach *i* drains to reach *i+1*). The network carries no
hydrological meaning -- it exists to exercise the coupling machinery. Checks
that runoff transfer is exact, drainage area accumulates correctly down the
chain, and mass is conserved through the network. This is the only script
here that needs nothing beyond a fresh checkout and a mizuRoute-enabled
build; run it after touching the coupling code.

## `test_mizuroute_provo_real_network.sh` -- the real Provo network

```
./test_mizuroute_provo_real_network.sh [work_dir]
```

Uses `hydrofabric_to_topology.py` to convert the **real** Provo hydrofabric
(`../test_ngen/domain_provo/settings/gage-10154200_subset.gpkg` -- the same
geopackage t-route reads directly) into a mizuRoute topology, then routes the
same SUMMA runoff through it. Checks the same mass-conservation properties as
the bundled test, on a real dendritic network with real reach lengths,
slopes, and drainage areas instead of a single synthetic chain.

`hydrofabric_to_topology.py` also cross-checks that the hydrofabric's
catchment areas agree with SUMMA's `attributes.nc` HRU areas (they're linked
only through each `cat-*.input` file's `attrib_file_HRU_order`, so this is a
real check, not a tautology).

## Comparing against t-route: `compare_to_troute.py`

`../test_ngen/provo_run.sh` runs the *same* Provo domain coupled inside ngen,
with t-route doing the routing instead of mizuRoute. Since both read the same
geopackage, mizuRoute's `segId` (as written by `hydrofabric_to_topology.py`)
and t-route's `feature_id` are literally the same hydrofabric catchment
number -- no id crosswalk needed.

To compare them:

1. `./test_mizuroute_provo_real_network.sh` (from this repo).
2. From the ngen repo root (needs ngen built and a python env with
   `nwm_routing` importable -- see `../test_ngen/readme.md`):
   `./extern/summa/summa/utils/test/test_ngen/provo_run.sh`
3. `python3 compare_to_troute.py <mizuroute .../output/*_timestep.nc> <ngen domain_provo/simulations dir with troute_output_*.nc> <mizuroute .../settings/topology.nc>`

**This is not expected to match closely, and that's expected, not a bug.**
mizuRoute is run here with the kinematic wave method (`methods = "3"` in
`test_mizuroute_provo_real_network.sh`), the same method the bundled test
already validates. t-route defaults to a Muskingum-Cunge-like scheme
(`compute_kernel: V02-structured` in `provo_routing.yaml`) that uses real
channel geometry -- width, Manning's n -- which the SUMMA-mizuRoute coupling
interface (`config.toml [hydrofabric]`, see `hydrofabric_to_topology.py`'s
docstring) does not currently pass through; only reach length, slope, and
drainage area are shared between the two runs. So `compare_to_troute.py`
is useful for confirming both are routing the same real network sensibly
(same reaches respond, comparable timing and volumes) -- it is not a
reference-quality validation of one against the other. Getting real
agreement would mean either running mizuRoute with its Muskingum-Cunge
method (`methods = "4"`) fed the hydrofabric's width/Manning's n -- which
means adding `varname_width`/`varname_man_n` wiring to the coupling's
`config.toml` parsing (`build/source/mizuroute/mizuroute_config.f90`) and to
`hydrofabric_to_topology.py` -- or configuring t-route's own kinematic wave
option to match mizuRoute instead.

## Files

- `make_test_topology.py` -- builds the synthetic network for
  `test_mizuroute_bundled.sh`.
- `hydrofabric_to_topology.py` -- converts a real NextGen hydrofabric
  geopackage to the same mizuRoute topology format, for
  `test_mizuroute_provo_real_network.sh`.
- `compare_to_troute.py` -- compares mizuRoute and t-route routed flow on a
  shared real network.

## Requirements

python3 with `netCDF4` and `numpy`. `hydrofabric_to_topology.py` reads the
geopackage directly via `sqlite3` (no `geopandas`/`fiona` needed).
