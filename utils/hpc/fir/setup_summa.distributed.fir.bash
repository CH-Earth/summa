#!/bin/bash

# Load modules and define environment variables
source "$HOME/models/summa/utils/hpc/fir/env_summa.fir.bash"

# ----------------------------------------------------------------------
# Define paths
# ----------------------------------------------------------------------

DATA_REL="projects/rrg-mclark/calibration/data/distributed/FH"
SOURCE_DIR="$HOME/$DATA_REL"
LOCAL_DIR="$SLURM_TMPDIR/data/distributed/FH"
DOMAIN_ARCHIVE="FH-domain.tar"

# ----------------------------------------------------------------------
# Stage SUMMA to node-local storage
# ----------------------------------------------------------------------

mkdir -p "$SLURM_TMPDIR/models/summa"

# Always refresh executables
rm -rf "$SLURM_TMPDIR/models/summa/bin"
cp -a "$HOME/models/summa/bin" \
      "$SLURM_TMPDIR/models/summa/"

# ----------------------------------------------------------------------
# Stage model inputs
# ----------------------------------------------------------------------

mkdir -p "$LOCAL_DIR"

if [ ! -d "$LOCAL_DIR/domain" ]; then

  echo "[$(hostname)] Staging distributed FH data..."

  cp "$SOURCE_DIR/$DOMAIN_ARCHIVE" \
     "$LOCAL_DIR/"

  tar --warning=no-unknown-keyword \
      -xf "$LOCAL_DIR/$DOMAIN_ARCHIVE" \
      -C "$LOCAL_DIR/"

  rm "$LOCAL_DIR/$DOMAIN_ARCHIVE"

  echo "[$(hostname)] Distributed FH data staged"

else

  echo "[$(hostname)] Distributed FH data already staged"

fi

# ----------------------------------------------------------------------
# Refresh experiment inputs
# ----------------------------------------------------------------------

for EXP_DIR in "$SOURCE_DIR"/exp*; do

  [ -d "$EXP_DIR" ] || continue

  EXP=$(basename "$EXP_DIR")
  SCRATCH_DIR="$LOCAL_DIR/$EXP"

  mkdir -p "$SCRATCH_DIR"

  echo "[$(hostname)] Refreshing common inputs for $EXP..."

  rm -rf "$SCRATCH_DIR/common_inputs"
  cp -a "$EXP_DIR/common_inputs" "$SCRATCH_DIR/"

  echo "[$(hostname)] Refreshing configuration files for $EXP..."

  rm -rf "$SCRATCH_DIR/settings"
  cp -a "$EXP_DIR/settings" "$SCRATCH_DIR/"

done

# ----------------------------------------------------------------------
# Stage observations
# ----------------------------------------------------------------------

OBS_SOURCE="$HOME/data/distributed/obs"
OBS_LOCAL="$SLURM_TMPDIR/data/distributed/obs"

mkdir -p "$OBS_LOCAL"

echo "[$(hostname)] Staging observations..."

cp -a "$OBS_SOURCE/CAN_07DD001_daily_flow_observations.nc" \
      "$OBS_LOCAL/"

echo "[$(hostname)] Observations staged"

# ----------------------------------------------------------------------
# Create runtime output directories
# ----------------------------------------------------------------------

for case_dir in "$LOCAL_DIR"/domain/*; do
  [ -d "$case_dir" ] || continue
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
echo "Local data directory: $LOCAL_DIR"
echo "Available CPUs:       $SLURM_CPUS_ON_NODE"
