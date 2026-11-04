# Multi-case calibration on Fir

Scripts for running the Century multi-basin calibration (PR #649) on the Alliance Fir cluster.
They assume SUMMA is cloned at `$HOME/models/summa`, SUNDIALS is installed at
`$HOME/local/sundials`, and the Century data is under `$HOME/data`. The data is not in this
repository; to exercise the same code path locally, use
`utils/test/test_calibration/multi_case_example/make_stub_century.bash`.

| File | Use |
| --- | --- |
| `build.fir.bash` | load modules and build SUMMA with MPI, SUNDIALS, mizuRoute and OpenMP |
| `env_summa.fir.bash` | runtime modules, and one GRU thread per core given to each MPI rank |
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

## Large domains

A trial's GRUs can run on `k` threads, and a case can span nodes with `nodes_per_case` in the manifest:

```bash
# 5 threads per trial, one case over 4 nodes: nodes_per_case = 4, cases_per_node = 1
sbatch --nodes=4 --ntasks=152 --ntasks-per-node=38 --cpus-per-task=5 summa_calib.fir.slurm <manifest>
```

Threads only help a generation if the workers still cover the population, since a generation takes
`ceil(population / workers)` trial times. Otherwise fewer threads and more trials at once is faster, so
one case of population 120 on one node is best at `k = 1`.

A trial cannot go faster than its largest GRU: the 3396-HRU Sagehen grid has 722 HRUs in GRU 1, which
bounds it at 4.7x. Measured there on an 8-core laptop: 1.5x at 2 threads, 2.0x at 3, 2.8x at 5, 3.0x at 9.
Results are bit-identical at any thread count.
