#!/bin/bash

# Load modules and define environment variables
source "$HOME/models/summa/build/cmake/env_summa.fir.bash"

# ----------------------------------------------------------------------
# Stage SUMMA to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/models/summa"

# Always refresh executables
rm -rf "$SLURM_TMPDIR/models/summa/bin"
cp -a "$HOME/models/summa/bin" \
      "$SLURM_TMPDIR/models/summa/"

# Always refresh calibration utilities
mkdir -p "$SLURM_TMPDIR/models/summa/utils/test"

rm -rf "$SLURM_TMPDIR/models/summa/utils/test/test_calibration"
cp -a "$HOME/models/summa/utils/test/test_calibration" \
      "$SLURM_TMPDIR/models/summa/utils/test/"

# ----------------------------------------------------------------------
# Stage model inputs
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/century/MM"

if [ ! -d "$SLURM_TMPDIR/data/century/MM/domain" ]; then

  echo "[$(hostname)] Staging Century data..."

  cp "$HOME/data/century/MM/MM-exp.tar" \
     "$SLURM_TMPDIR/data/century/MM/"

  tar --warning=no-unknown-keyword \
      -xf "$SLURM_TMPDIR/data/century/MM/MM-exp.tar" \
      -C "$SLURM_TMPDIR/data/century/MM/"

  rm "$SLURM_TMPDIR/data/century/MM/MM-exp.tar"

  echo "[$(hostname)] Century data staged"

else

  echo "[$(hostname)] Century data already staged"

fi

# ----------------------------------------------------------------------
# Stage observations
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/data/camels-spat/observations"

if [ ! -d "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily" ]; then

  echo "[$(hostname)] Staging observations..."

  cp "$HOME/data/camels-spat/observations/obs-daily.tar" \
     "$SLURM_TMPDIR/data/camels-spat/observations/"

  tar --warning=no-unknown-keyword \
      -xf "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar" \
      -C "$SLURM_TMPDIR/data/camels-spat/observations/"

  rm "$SLURM_TMPDIR/data/camels-spat/observations/obs-daily.tar"

  echo "[$(hostname)] Observations staged"

else

  echo "[$(hostname)] Observations already staged"

fi

# ----------------------------------------------------------------------
# Create runtime output directories
# ----------------------------------------------------------------------

for case_dir in "$SLURM_TMPDIR"/data/century/MM/domain/*; do
  mkdir -p "$case_dir/work"
done

# ----------------------------------------------------------------------
# Move to node-local working directory
# ----------------------------------------------------------------------

cd "$SLURM_TMPDIR/models/summa"

echo
echo "SUMMA runtime environment ready"
echo "Host:                 $(hostname)"
echo "Local work directory: $SLURM_TMPDIR/models/summa"
echo "Local data directory: $SLURM_TMPDIR/data"
echo "Available CPUs:       $SLURM_CPUS_ON_NODE"
