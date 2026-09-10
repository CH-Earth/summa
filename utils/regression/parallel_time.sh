#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

OUTPUT_DIR="$HOME/data/great-slave-lake/Athabasca-river/results/SUMMA"
ORIGINAL_DIR="$OUTPUT_DIR/original"

NP_LIST="2 4 8 12 16"

TMP_DIR="$OUTPUT_DIR/wallclock_tmp"
mkdir -p "$TMP_DIR"


# ==============================================================================
# Helper: total wall-clock time in one output file
# ==============================================================================

get_runtime () {

  local file=$1
  local tmp="$TMP_DIR/runtime.nc"

  ncap2 -O \
    -s 'runtime=wallClockTime.total()' \
    "$file" "$tmp"

  ncks -H -C -s '%g' -v runtime "$tmp"
}


# ==============================================================================
# Helper: summarize one set of MPI experiments
# ==============================================================================

summarize_scaling () {

  local output_dir=$1
  local label=$2

  local np1_file="$output_dir/run_1_dev_np1_timestep.nc"
  local np1_time

  np1_time=$(get_runtime "$np1_file")

  echo
  echo "$label"
  echo "========================================================================="
  printf "%-12s %10s %10s %12s\n" "Experiment" "Time (s)" "Speedup" "Efficiency"
  echo "-------------------------------------------------------------------------"

  printf "%-12s %10.3f %10.2f %11.1f%%\n" \
    "np1" "$np1_time" 1.0 100.0

  for np in $NP_LIST; do

    max_time=0

    for file in "$output_dir"/run_1_dev_np${np}_G*-*_timestep.nc; do

      runtime=$(get_runtime "$file")

      if (( $(echo "$runtime > $max_time" | bc -l) )); then
        max_time=$runtime
      fi

    done

    speedup=$(echo "$np1_time / $max_time" | bc -l)
    efficiency=$(echo "100.0 * $speedup / $np" | bc -l)

    printf "%-12s %10.3f %10.2f %11.1f%%\n" \
      "np${np}" "$max_time" "$speedup" "$efficiency"

  done
}


# ==============================================================================
# Timing summaries
# ==============================================================================

summarize_scaling "$OUTPUT_DIR"   "New implementation"
summarize_scaling "$ORIGINAL_DIR" "Original MPI implementation"


# ==============================================================================
# Cleanup
# ==============================================================================

rm -rf "$TMP_DIR"
