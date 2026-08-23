#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

OUTPUT_DIR="$HOME/data/great-slave-lake/Athabasca-river/results/SUMMA"

STABLE_FILE="$OUTPUT_DIR/run_1_stable_serial_timestep.nc"
DEV_FILE="$OUTPUT_DIR/run_1_dev_serial_timestep.nc"
MPI1_FILE="$OUTPUT_DIR/run_1_dev_np1_timestep.nc"

# HRU-level variables
HRU_VARIABLES=(
  scalarCanopyTemp
  scalarSurfaceTemp
  scalarRootZoneTemp
  scalarCanopyWat
  scalarSWE
  scalarTotalSoilWat
)

# GRU-level variables
GRU_VARIABLES=(
  averageRoutedRunoff
)

hru_vars=$(IFS=,; echo "${HRU_VARIABLES[*]}")
gru_vars=$(IFS=,; echo "${GRU_VARIABLES[*]}")

TMP_DIR="$OUTPUT_DIR/regression_tmp"
mkdir -p "$TMP_DIR"

# summary arrays
SUMMARY_LABEL=()
SUMMARY_DIFF=()


# ==============================================================================
# Helper: compare variables in two files
# ==============================================================================

compare_files () {

  file1=$1
  file2=$2
  label=$3
  shift 3

  variables=("$@")
  var_list=$(IFS=,; echo "${variables[*]}")

  diff_file="$TMP_DIR/diff.nc"
  summary_file="$TMP_DIR/summary.nc"

  LAST_MAX_DIFF=0

  echo
  echo "$label"
  echo "----------------------------------------------------------------"

  ncdiff -O -v "$var_list" \
    "$file1" "$file2" \
    "$diff_file"

  for var in "${variables[@]}"; do

    ncap2 -O \
      -s "max_diff=max(abs(${var}))" \
      "$diff_file" "$summary_file"

    max_diff=$(ncks -H -C -s '%g' -v max_diff "$summary_file")

    printf "%-30s max_abs_diff = %s\n" "$var" "$max_diff"

    if (( $(echo "$max_diff > $LAST_MAX_DIFF" | bc -l) )); then
      LAST_MAX_DIFF=$max_diff
    fi

  done
}


# ==============================================================================
# 1. Stable serial vs development serial
# ==============================================================================

compare_files \
  "$STABLE_FILE" "$DEV_FILE" \
  "stable serial vs development serial" \
  "${HRU_VARIABLES[@]}" "${GRU_VARIABLES[@]}"

SUMMARY_LABEL+=("stable serial vs dev serial")
SUMMARY_DIFF+=("$LAST_MAX_DIFF")


# ==============================================================================
# 2. Development serial vs MPI np=1
# ==============================================================================

compare_files \
  "$DEV_FILE" "$MPI1_FILE" \
  "development serial vs MPI np=1" \
  "${HRU_VARIABLES[@]}" "${GRU_VARIABLES[@]}"

SUMMARY_LABEL+=("dev serial vs MPI np=1")
SUMMARY_DIFF+=("$LAST_MAX_DIFF")


# ==============================================================================
# 3. Development serial vs decomposed MPI runs
# ==============================================================================

for np in 2 4 8 12 16; do

  np_max_diff=0

  echo
  echo
  echo "================================================================"
  echo "development serial vs MPI np=$np"
  echo "================================================================"

  for mpi_file in "$OUTPUT_DIR"/run_1_dev_np${np}_G*-*_timestep.nc; do

    # --------------------------------------------------------------------------
    # determine the HRU and GRU ranges represented in this processor file
    # --------------------------------------------------------------------------

    range_file="$TMP_DIR/range.nc"

    ncap2 -O \
      -s 'hru_min=min(hru); hru_max=max(hru); gru_min=min(gru); gru_max=max(gru)' \
      "$mpi_file" "$range_file"

    hru_min=$(ncks -H -C -s '%d' -v hru_min "$range_file")
    hru_max=$(ncks -H -C -s '%d' -v hru_max "$range_file")

    gru_min=$(ncks -H -C -s '%d' -v gru_min "$range_file")
    gru_max=$(ncks -H -C -s '%d' -v gru_max "$range_file")

    echo
    echo "$(basename "$mpi_file")"
    echo "  HRUs: $hru_min-$hru_max"
    echo "  GRUs: $gru_min-$gru_max"

    # --------------------------------------------------------------------------
    # HRU variables
    # --------------------------------------------------------------------------

    serial_subset="$TMP_DIR/serial_hru.nc"

    ncks -F -O \
      -v "$hru_vars" \
      -d hru,"$hru_min","$hru_max" \
      "$DEV_FILE" "$serial_subset"

    compare_files \
      "$serial_subset" "$mpi_file" \
      "  HRU variables" \
      "${HRU_VARIABLES[@]}"

    if (( $(echo "$LAST_MAX_DIFF > $np_max_diff" | bc -l) )); then
      np_max_diff=$LAST_MAX_DIFF
    fi

    # --------------------------------------------------------------------------
    # GRU variables
    # --------------------------------------------------------------------------

    serial_subset="$TMP_DIR/serial_gru.nc"

    ncks -F -O \
      -v "$gru_vars" \
      -d gru,"$gru_min","$gru_max" \
      "$DEV_FILE" "$serial_subset"

    compare_files \
      "$serial_subset" "$mpi_file" \
      "  GRU variables" \
      "${GRU_VARIABLES[@]}"

    if (( $(echo "$LAST_MAX_DIFF > $np_max_diff" | bc -l) )); then
      np_max_diff=$LAST_MAX_DIFF
    fi

  done

  SUMMARY_LABEL+=("dev serial vs MPI np=$np")
  SUMMARY_DIFF+=("$np_max_diff")

done


# ==============================================================================
# 4. Regression summary
# ==============================================================================

echo
echo
echo "Regression summary"
echo "======================================================================"
printf "%-32s %16s %8s\n" "Comparison" "Max abs diff" "Result"
echo "----------------------------------------------------------------------"

for ((i=0; i<${#SUMMARY_LABEL[@]}; i++)); do

  diff=${SUMMARY_DIFF[$i]}

  if (( $(echo "$diff == 0" | bc -l) )); then
    result="PASS"
  else
    result="FAIL"
  fi

  printf "%-32s %16s %8s\n" \
    "${SUMMARY_LABEL[$i]}" "$diff" "$result"

done


# ==============================================================================
# Cleanup
# ==============================================================================

rm -rf "$TMP_DIR"
