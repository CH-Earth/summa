#!/bin/bash
set -euo pipefail

# load modules and define environment variables
source "$HOME/models/summa/build/cmake/env_summa.fir.bash"

# ----------------------------------------------------------------------
# Stage SUMMA executable to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/models/summa"

echo "[$(hostname)] Refreshing SUMMA executable..."

rm -rf "$SLURM_TMPDIR/models/summa/bin"

cp -a "$HOME/models/summa/bin" \
      "$SLURM_TMPDIR/models/summa/"

echo "[$(hostname)] SUMMA executable refreshed"

# ----------------------------------------------------------------------
# Stage recession domains to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/recession"

if [ ! -d "$SLURM_TMPDIR/data/recession/domain" ]; then

  echo "[$(hostname)] Copying recession domain archive..."

  cp "$HOME/data/recession/domain.tar" \
     "$SLURM_TMPDIR/data/recession/"

  echo "[$(hostname)] Extracting recession domains..."

  tar --warning=no-unknown-keyword \
      -xf "$SLURM_TMPDIR/data/recession/domain.tar" \
      -C "$SLURM_TMPDIR/data/recession/"

  rm "$SLURM_TMPDIR/data/recession/domain.tar"

  echo "[$(hostname)] Recession domains extracted"

else

  echo "[$(hostname)] Recession domains already staged"

fi

# ----------------------------------------------------------------------
# Refresh shared recession experiment inputs
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/recession/exp01"

echo "[$(hostname)] Refreshing common inputs..."

rm -rf "$SLURM_TMPDIR/data/recession/exp01/common_inputs"

cp -a "$HOME/data/recession/exp01/common_inputs" \
      "$SLURM_TMPDIR/data/recession/exp01/"

echo "[$(hostname)] Refreshing configuration files..."

rm -rf "$SLURM_TMPDIR/data/recession/exp01/settings"

cp -a "$HOME/data/recession/exp01/settings" \
      "$SLURM_TMPDIR/data/recession/exp01/"

# ----------------------------------------------------------------------
# Stage observations to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/camels-spat/observations"

if [ ! -d "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily" ]; then

  echo "[$(hostname)] Copying observations archive..."

  cp "$HOME/data/camels-spat/observations/obs-daily.tar" \
     "$SLURM_TMPDIR/data/camels-spat/observations/"

  echo "[$(hostname)] Extracting observations..."

  tar --warning=no-unknown-keyword \
      -xf "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar" \
      -C "$SLURM_TMPDIR/data/camels-spat/observations/"

  rm "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar"

  echo "[$(hostname)] Observations extracted"

else

  echo "[$(hostname)] Observations already staged"

fi

# ----------------------------------------------------------------------
# Prepare runtime output directories
# ----------------------------------------------------------------------

for case_dir in "$SLURM_TMPDIR"/data/recession/domain/*; do

  [ -d "$case_dir" ] || continue

  mkdir -p "$case_dir/work"

done

# ----------------------------------------------------------------------
# Move to node-local working directory
# ----------------------------------------------------------------------

cd "$SLURM_TMPDIR/models/summa"

echo
echo "SUMMA recession runtime environment ready"
echo "Host:                  $(hostname)"
echo "Local model directory: $SLURM_TMPDIR/models/summa"
echo "Local data directory:  $SLURM_TMPDIR/data"
echo "Available CPUs:         $SLURM_CPUS_ON_NODE"
echo
