#!/usr/bin/env bash
set -euo pipefail

# Command-line arguments, with environment-variable fallbacks
RUN_LABEL="${1:-${RUN_LABEL:-}}"
SUMMA_DATA_ROOT="${2:-${SUMMA_DATA_ROOT:-}}"

: "${SLURM_JOB_ID:?Must run inside a Slurm allocation}"
: "${RUN_LABEL:?Usage: $0 <run_label> <data_root>}"
: "${SUMMA_DATA_ROOT:?Shared dataset root must be specified}"

export RUN_LABEL SUMMA_DATA_ROOT

# Initialize node-local storage
export SLURM_TMPDIR="${SLURM_TMPDIR:-/tmp/${USER}_${SLURM_JOB_ID}}"

mkdir -p "$SLURM_TMPDIR"
chmod 700 "$SLURM_TMPDIR"

echo "[$(hostname)] SLURM_TMPDIR=$SLURM_TMPDIR"

REPO="$HOME/models/summa"
SHARED_MM="$SUMMA_DATA_ROOT/century/MM"
SHARED_OBS="$SUMMA_DATA_ROOT/camels-spat/observations"
LOCAL_MM="$SLURM_TMPDIR/data/century/MM"
LOCAL_OBS="$SLURM_TMPDIR/data/camels-spat/observations"
LOCAL_REPO="$SLURM_TMPDIR/models/summa"

DOMAIN_TAR="$SHARED_MM/MM-domains.tar"
OBS_TAR="$SHARED_OBS/obs-daily.tar"
SHARED_EXP="$SHARED_MM/$RUN_LABEL"

for item in "$DOMAIN_TAR" "$OBS_TAR" \
            "$SHARED_EXP/common_inputs" "$SHARED_EXP/settings" \
            "$SHARED_EXP/settings/manifest_century.toml" \
            "$REPO/bin/summa_sundials_mizuroute_opt.exe"; do
  if [[ ! -e "$item" ]]; then
    echo "[$(hostname)] ERROR: Missing required input: $item" >&2
    exit 1
  fi
done

mkdir -p "$SLURM_TMPDIR"
chmod 700 "$SLURM_TMPDIR"
mkdir -p "$LOCAL_REPO" "$LOCAL_MM" "$LOCAL_OBS"

# Refresh the executable on each allocated compute node.
rm -rf "$LOCAL_REPO/bin"
cp -a "$REPO/bin" "$LOCAL_REPO/"

# Avoid doubling peak /tmp usage: extract directly from the shared archive.
echo "[$(hostname)] Extracting domains from $DOMAIN_TAR"
tar --warning=no-unknown-keyword -xf "$DOMAIN_TAR" -C "$LOCAL_MM"
if [[ ! -d "$LOCAL_MM/domain" ]]; then
  echo "[$(hostname)] ERROR: Domain archive must extract domain/<case>/..." >&2
  exit 1
fi

# Only stage the chosen experiment (the Fir script staged every exp*).
echo "[$(hostname)] Staging experiment: $RUN_LABEL"
mkdir -p "$LOCAL_MM/$RUN_LABEL"
cp -a "$SHARED_EXP/common_inputs" "$LOCAL_MM/$RUN_LABEL/"
cp -a "$SHARED_EXP/settings" "$LOCAL_MM/$RUN_LABEL/"

# Stage the observations used by the Century workflow.
echo "[$(hostname)] Extracting observations from $OBS_TAR"
tar --warning=no-unknown-keyword -xf "$OBS_TAR" -C "$LOCAL_OBS"
if [[ ! -d "$LOCAL_OBS/obs-daily" ]]; then
  echo "[$(hostname)] ERROR: Observation archive must extract obs-daily/" >&2
  exit 1
fi

shopt -s nullglob
for case_dir in "$LOCAL_MM"/domain/*; do
  [[ -d "$case_dir" ]] && mkdir -p "$case_dir/work"
done

echo "[$(hostname)] Local data staged at $SLURM_TMPDIR"
df -h "$SLURM_TMPDIR"
