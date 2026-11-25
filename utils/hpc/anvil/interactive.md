## Steps to run an interactive run

### Load SUMMA environment
```bash
cd ~/models/summa
source utils/hpc/anvil/env_summa.anvil.bash
```

### Stage SUMMA and calibration data
```bash
export DATA_ROOT=/anvil/projects/x-ees250064/summa-calibration/calibration/data
bash utils/hpc/anvil/setup_summa.century.anvil.bash exp50 $DATA_ROOT
```
### Define the local path in your interactive shell
```bash
export SLURM_TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
cd "$SLURM_TMPDIR/models/summa"
```

### Create a directory for persistent output
```bash
export OUT="$DATA_ROOT/../test-results/exp50_smoke_${SLURM_JOB_ID}"
mkdir -p "$OUT"
```

### If necessary, modify control files
```bash
mkdir -p ~/summa-tests/exp50
cp $DATA_ROOT/century/MM/exp50/settings/manifest_century.toml ~/summa-tests/exp50/manifest_test.toml
cp $DATA_ROOT/century/MM/exp50/settings/summa_config_template.toml ~/summa-tests/exp50/summa_config_test.toml
```
(and modify as appropriate)

### Run SUMMA (using the manifest in the $HOME directory)
```bash
srun -n 4 \
  ./bin/summa_sundials_mizuroute_opt.exe \
  --manifest "$HOME/summa-tests/exp50/manifest_test.toml" \
  --home "$SLURM_TMPDIR"
```
