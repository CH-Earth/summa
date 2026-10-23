# test_calibration

Tests and worked examples for SUMMA's parameter calibration (`summa[_sundials]_opt.exe`).

## Requirements

Calibration needs an executable built with **both** `-DUSE_MPI=ON` and `-DUSE_MIZUROUTE=ON`.

MPI is what builds the calibration driver in the first place, and it is how trials are run
concurrently: one rank coordinates, the others evaluate parameter sets.

mizuRoute is needed because the objective function compares *routed streamflow* against gauge
observations, and routed streamflow is produced by mizuRoute. In a build without it the
simulated flow series is never filled, so there is nothing to score. This is a real constraint
of the current design, not a property of these tests.

## `test_calibration_bow.sh`

A short, self-contained calibration on the Bow River at Banff (CAN_05BB001). Everything it needs
is already in the repository, under `utils/test/test_mizuroute/bow_real_data/`: lumped SUMMA
inputs and ten years of forcing, the mizuRoute topology and lumped-to-HRU remapping, and daily
streamflow observations.

```bash
./test_calibration_bow.sh [n_samples] [n_ranks]     # defaults: 4 samples, 2 ranks
```

It writes a TOML configuration pointing at that data, runs a DDS search over five parameters for
one simulated year after a year of spinup, and checks that finite, plausible objective values were
recorded for the trials.

The period and the sample count are deliberately tiny so the test finishes in minutes. It
exercises the machinery; it does **not** produce a calibrated parameter set. A real calibration
uses thousands of samples over several years — see `n_samples` in the generated config.

## Choosing the parameters to calibrate

List them in `[calibration.parameters]`, one entry per parameter. Everything in an entry is optional:

```toml
[calibration.parameters]
k_soil             = { lower = 1.0e-7, upper = 1.0e-4, transformation = "log" }
theta_sat          = {}                       # native bounds from the parameter files, no transformation
aquiferScaleFactor = { lower = 0.15 }         # own lower bound, native upper bound
```

A bound that is left out is taken from the SUMMA parameter files (`localParamInfo.txt`,
`basinParamInfo.txt`), or from the code defaults if they do not list it. The transformation
defaults to `none`. A bound set outside the native range is allowed but logs a warning.

The older form still works, alone or together with the table: `param_list = ["k_soil", ...]` in
`[calibration]`, with transformations given in `[calibration.parameter_transformations]`
(`k_soil = "log"`).

## `multi_case_example/` -- multi-case calibration

`--manifest <file>` calibrates many basins in one job, from a manifest listing the cases and a
per-case configuration template that the manifest expands.

The test dataset for this is built by `make_stub_century.bash`, which assembles the directory
layout the template expects from the one domain bundled in this repository, the Bow at Banff.
Files are symlinked, so a three-case dataset costs about 128 KB:

```bash
cd multi_case_example
./make_stub_century.bash                  # defaults: stub_century/, 3 cases, 1 at a time
mpirun -np 2 ../../../../bin/summa_sundials_mizuroute_opt.exe \
       --manifest stub_century/stub_manifest.toml
```

Arguments are `[root] [n_cases] [cases_per_node]`. Each case group needs at least two ranks --
one coordinates the search, the rest evaluate samples -- and the groups have to divide the node's
ranks evenly. If the rank count does not fit, the driver warns and reduces `cases_per_node` rather
than refusing to start, so `-np` equal to `2 * cases_per_node` gets the cases you asked for and
anything else still runs. The script prints a matching command when it finishes.

Every case is the same basin under a different name, so the calibrated values are not meaningful.
What it exercises is the part with no other coverage: reading the manifest, expanding the template
per case, distributing cases across ranks, and writing per-case output. Two real bugs in the
multi-case path were found this way, both fixed -- calibration structures were not released
between cases, so the second case failed to allocate, and the template's `work_path` lacked a
trailing slash, so output landed beside the work directory rather than in it.

To run it on a larger dataset, point `home_path` in the template at that dataset and list its
cases in the manifest; the layout is the same.

| File | What it is |
| --- | --- |
| `make_stub_century.bash` | builds the stub dataset, and the manifest and template that drive it |
| `manifest_century.toml` | the full manifest, listing ~115 basins; the stub manifest is derived from it |
| `summa_config_template.toml` | the per-case configuration template the manifest expands |
| `summa_config_CAN_05BB001.toml` | a single filled-in case, useful for seeing what the template produces |
| `century_cases.txt` | the basin list |
| `setup_summa_cases.bash` | builds per-case input directories from a larger dataset, if you have one |
| `check_time.bash` | reports each case's forcing time range, to check a calibration period is covered |

`setup_summa_cases.bash` and `check_time.bash` hard-code a data root (`$HOME/data/century/...`)
and an experiment-specific layout, so they need editing before they run anywhere else. They are
kept because they document how the per-case directories were originally produced;
`make_stub_century.bash` is what builds a dataset you can run today.
