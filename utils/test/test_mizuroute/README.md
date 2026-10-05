# SUMMA-mizuRoute coupling tests

Tests for the built-in SUMMA-mizuRoute coupling (`-DUSE_MIZUROUTE=ON`, see
[docs/index.md](../../../docs/index.md)). There are three, unrelated except
for all exercising the same coupling code -- pick the one that matches what
you're checking:

| Test case | What it is | Run with | Output goes to |
|---|---|---|---|
| **toy problem** | Bundled Provo domain + a made-up, single-chain river network. No real-world meaning; exists purely to exercise the coupling machinery. | `./test_mizuroute_bundled.sh` | `toy_prob/` (generated, gitignored) |
| **Provo, real network** | Same bundled Provo domain, but routed through its *real* hydrofabric network, so it can be compared against the same domain routed with t-route/ngen. | `./test_mizuroute_provo_real_network.sh`, then `compare_to_troute.py` | `provo_real_network/` (generated, gitignored) |
| **Bow real data** | A real, independent SUMMA-mizuRoute case study (Bow River at Banff): real forcing, real network + remapping, observed streamflow, KGE evaluation. Not a pass/fail check against anything above. | see [`bow_real_data/README.md`](bow_real_data/README.md), then `check_objective_func.R` | `bow_real_data/work/` (generated, gitignored) |

Only the toy problem is expected to run unmodified right after a fresh
checkout -- it needs nothing beyond a mizuRoute-enabled build. `bow_real_data`
is a tracked input dataset (real data, checked into git), not a generated
directory; the other two work dirs only appear once you run their script.

## `test_mizuroute_bundled.sh` -- toy problem: synthetic network, fully self-contained

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

**A close match is not expected, and that's expected, not a bug** -- but it
should be a *reasonably* close one (median NSE ~0.98 on this domain), not
wildly off. If you see something like a systematic ~3600x scale difference,
or several reaches pinned at exactly zero t-route flow for the whole run,
that's not this -- see the two coupling details below.

- **Units.** ngen (`include/core/Layer.hpp`) and t-route both assume a
  module's main output variable is a runoff rate in m/h, the convention
  CFE, LSTM and dHBV2 follow. SUMMA's BMI therefore also exposes
  `land_surface_water__runoff_volume_flux_mh` (the same runoff as
  `land_surface_water__runoff_volume_flux`, in `m h-1` instead of `m s-1`),
  and the realization configs use it as `main_output_variable`. Using the
  m/s variable instead makes routed flow ~3600x too small.
- **Confluences.** At a nexus fed by more than one catchment (any
  tributary junction), ngen sums the contributing catchments' flows into
  the single `nex-*_output.csv` value, so t-route can't tell which
  flowpath each part belongs to, and several reaches get zero flow. The
  provo routing config (`provo_routing.yaml`) therefore reads the
  per-catchment `cat-*.csv` files instead
  (`qlat_file_pattern_filter: "cat-*"`, `qlat_file_value_col:
  land_surface_water__runoff_volume_flux_mh`). This needs a t-route that
  supports `cat-*` input (CIROH-UA/t-route `ngiab`, PR #29).

With both handled, what's left is genuine routing-scheme difference: mizuRoute
is run here with the kinematic wave method (`methods = "3"` in
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

## `check_objective_func.R` -- validating the objective-function code, on Bow real data

Independently recomputes KGE/NSE/RMSE/MAE in R from `bow_real_data`'s routed
streamflow and observations, and cross-checks against what SUMMA itself
computed: run `bow_real_data` with `write_aligned = true` under `[objective]`
in its TOML (see [`bow_real_data/README.md`](bow_real_data/README.md)), and
this script also reads back the `objective`/`eval_qobs`/`eval_qsim` values
`write_evaluation()` (`build/source/objfunc/write_evaluation.f90`) wrote into
the same output file, and plots both alignments together. Needs
`bow_real_data` to have been run first; edit the evaluation period at the top
if you change the TOML's.

## Private-data mizuRoute validation: `compare_mizuroute.R`, `compare_coupled_sequential_mizuroute.R`

Two more R scripts doing coupled-vs-standalone mizuRoute comparisons like
`test_mizuroute_coupling.sh`/`../test_regression/`, plotted rather than
diffed, against private data (edit the hardcoded paths at the top before
running -- neither runs against anything bundled in this repo):

- `compare_mizuroute.R` -- three-way comparison (KGE/NSE, one plot) between
  standalone mizuRoute reading GRU-remapped runoff, standalone mizuRoute
  reading HRU runoff directly, and the SUMMA-mizuRoute coupled run, all on
  the same reach.
- `compare_coupled_sequential_mizuroute.R` -- coupled vs. sequential
  (standalone mizuRoute fed SUMMA's output) on the GSL/Athabasca domain, same
  idea as `test_mizuroute_coupling.sh` but as a KGE/NSE + plot instead of a
  reach-by-reach diff.

Both need a standalone mizuRoute executable, separate from SUMMA's own
build -- `make_mizuRoute.sh` builds one (macOS/Homebrew example; see the
[mizuRoute build docs](https://mizuroute.readthedocs.io/en/main/users_guide/Build_model.html)
for other platforms) and copies it to `bin/route_runoff`.
`test_mizuroute_coupling.sh` in `../test_regression/` needs the same thing.

## Files

- `make_test_topology.py` -- builds the synthetic network for
  `test_mizuroute_bundled.sh` (toy problem).
- `hydrofabric_to_topology.py` -- converts a real NextGen hydrofabric
  geopackage to the same mizuRoute topology format, for
  `test_mizuroute_provo_real_network.sh` (Provo, real network).
- `compare_to_troute.py` -- compares mizuRoute and t-route routed flow on a
  shared real network (Provo, real network).
- `bow_real_data/` -- tracked input data and settings for the Bow real data
  case study; see its own [README](bow_real_data/README.md).
- `check_objective_func.R`, `compare_mizuroute.R`,
  `compare_coupled_sequential_mizuroute.R` -- R comparison/validation
  scripts, described above.
- `make_mizuRoute.sh` -- builds the standalone `route_runoff` executable
  needed by `compare_mizuroute.R`, `compare_coupled_sequential_mizuroute.R`,
  and `../test_regression/test_mizuroute_coupling.sh`.

## Requirements

python3 with `netCDF4` and `numpy`. `hydrofabric_to_topology.py` reads the
geopackage directly via `sqlite3` (no `geopandas`/`fiona` needed). The R
scripts need `ncdf4` and `hydroGOF`, and source
`../test_regression/plot_utils.R`.
