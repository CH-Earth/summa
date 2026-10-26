# Multi-case calibration on Fir

Scripts for running the Century multi-basin calibration (PR #649) on the Alliance Fir cluster.
They assume SUMMA is cloned at `$HOME/models/summa`, SUNDIALS is installed at
`$HOME/local/sundials`, and the Century data is under `$HOME/data`. The data is not in this
repository; to exercise the same code path locally, use
`utils/test/test_calibration/multi_case_example/make_stub_century.bash`.

| File | Use |
| --- | --- |
| `build.fir.bash` | load modules and build SUMMA with MPI, SUNDIALS and mizuRoute |
| `env_summa.fir.bash` | runtime modules, and one compute thread per MPI rank |
| `setup_summa.century.fir.bash` | stage `bin/`, the Century data and observations to node-local `$SLURM_TMPDIR` |
| `summa_calib.fir.slurm` | stage once per node, then run the calibration from a manifest |

```bash
./build.fir.bash
sbatch summa_calib.fir.slurm                                    # 1 node, default manifest
sbatch --nodes=4 --ntasks=768 summa_calib.fir.slurm <manifest>  # 4 nodes
```

Data layout expected under `$HOME/data`, archived so staging is one copy and one extract:

```
century/MM/MM-exp.tar            -> century/MM/domain/<case>/...
century/MM/exp*/common_inputs/
century/MM/exp*/settings/        manifest and summa_config_template.toml for the run
camels-spat/observations/obs-daily.tar
```

The manifest and template in `settings/` follow
`utils/test/test_calibration/multi_case_example/`, with `basin_dir` and the common settings
pointed at the staged layout. `--home $SLURM_TMPDIR` roots every path at the node-local copy,
and `--persistent-output` copies each finished basin's calibration file to
`$SCRATCH/calibration/<job id>`.

Scaling is set by the node count and the manifest's `cases_per_node`; each case group gets
`192 / cases_per_node` ranks, one coordinating and the rest evaluating parameter sets. For
112 basins at 1000 evaluations each, 4 cases per node took 5.2 h on one node and 0.84 h on eight.
