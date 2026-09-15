#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

DATA_DIR="$HOME/data/great-slave-lake/Athabasca-river"

MIZU_FILE="$DATA_DIR/results/mizuRoute/run_1_standalone.h.2014-09-01-03600.nc"
SUMMA_FILE="$DATA_DIR/results/SUMMA/run_1_coupled_timestep.nc"

TMP_DIR="$DATA_DIR/results/mizuRoute/regression_tmp"

MIZU_TEST="$TMP_DIR/mizuroute.nc"
SUMMA_TEST="$TMP_DIR/summa_mizuroute.nc"
DIFF_FILE="$TMP_DIR/diff.nc"
SUMMARY_FILE="$TMP_DIR/summary.nc"


# ==============================================================================
# Setup
# ==============================================================================

echo
echo "======================================================================"
echo "mizuRoute coupling regression test"
echo "======================================================================"
echo
echo "Standalone mizuRoute:"
echo "  $MIZU_FILE"
echo
echo "Coupled SUMMA-mizuRoute:"
echo "  $SUMMA_FILE"
echo
echo "Temporary directory:"
echo "  $TMP_DIR"
echo

mkdir -p "$TMP_DIR"


# ==============================================================================
# Check input files
# ==============================================================================

echo "Checking input files..."

if [[ ! -f "$MIZU_FILE" ]]; then
  echo "ERROR: standalone mizuRoute file not found"
  echo "  $MIZU_FILE"
  exit 1
fi

if [[ ! -f "$SUMMA_FILE" ]]; then
  echo "ERROR: coupled SUMMA file not found"
  echo "  $SUMMA_FILE"
  exit 1
fi

echo "  input files found"
echo


# ==============================================================================
# Prepare standalone mizuRoute output
# ==============================================================================

echo "Preparing standalone mizuRoute output..."

ncks -O \
  -v time,reachID,KWroutedRunoff \
  "$MIZU_FILE" \
  "$MIZU_TEST"

echo "  wrote:"
echo "  $MIZU_TEST"
echo


# ==============================================================================
# Prepare coupled SUMMA-mizuRoute output
# ==============================================================================

echo "Preparing coupled SUMMA-mizuRoute output..."

echo "  extracting q_reach..."

ncks -O \
  -v time,seg,q_reach \
  -d method,0 \
  "$SUMMA_FILE" \
  "$SUMMA_TEST"

echo "  removing singleton method dimension..."

ncwa -O \
  -a method \
  "$SUMMA_TEST" \
  "$SUMMA_TEST"

echo "  converting q_reach to standalone mizuRoute precision..."

ncap2 -O \
  -s 'KWroutedRunoff=float(q_reach)' \
  "$SUMMA_TEST" \
  "$SUMMA_TEST"

echo "  retaining comparison variables..."

ncks -O \
  -v time,seg,KWroutedRunoff \
  "$SUMMA_TEST" \
  "$SUMMA_TEST"

echo "  renaming seg to reachID..."

ncrename -O \
  -v seg,reachID \
  "$SUMMA_TEST"

echo "  wrote:"
echo "  $SUMMA_TEST"
echo


# ==============================================================================
# Compare reach IDs
# ==============================================================================

echo "Comparing reach IDs..."

MIZU_IDS="$TMP_DIR/reach_mizu.txt"
SUMMA_IDS="$TMP_DIR/reach_summa.txt"

ncks -H -C -s '%d\n' -v reachID "$MIZU_TEST"  > "$MIZU_IDS"
ncks -H -C -s '%d\n' -v reachID "$SUMMA_TEST" > "$SUMMA_IDS"

# exact same set of reach IDs
sort -n "$MIZU_IDS"  > "$TMP_DIR/reach_mizu_sorted.txt"
sort -n "$SUMMA_IDS" > "$TMP_DIR/reach_summa_sorted.txt"

if cmp -s "$TMP_DIR/reach_mizu_sorted.txt" \
          "$TMP_DIR/reach_summa_sorted.txt"; then
  echo "  PASS: reach ID sets match exactly"
else
  echo "  FAIL: reach ID sets differ"
  diff "$TMP_DIR/reach_mizu_sorted.txt" \
       "$TMP_DIR/reach_summa_sorted.txt" || true
  exit 1
fi

echo

# ==============================================================================
# Compare routed streamflow reach by reach
# ==============================================================================

n_reaches=0
n_pass=0
n_fail=0
max_diff_reach=""

echo "Comparing routed streamflow by reach..."

max_diff=0
flow_result="PASS"

for reach_id in $(cat "$MIZU_IDS"); do

  i_mizu=$(awk -v id="$reach_id" '$1==id {print NR-1}' "$MIZU_IDS")
  i_summa=$(awk -v id="$reach_id" '$1==id {print NR-1}' "$SUMMA_IDS")

  mizu_reach="$TMP_DIR/mizu_reach.nc"
  summa_reach="$TMP_DIR/summa_reach.nc"
  diff_reach="$TMP_DIR/diff_reach.nc"
  summary_reach="$TMP_DIR/summary_reach.nc"

  # ----- standalone mizuRoute: KWroutedRunoff -----

  ncks -O \
    -v KWroutedRunoff \
    -d seg,"$i_mizu","$i_mizu" \
    "$MIZU_FILE" \
    "$mizu_reach"

  # ----- coupled SUMMA-mizuRoute: q_reach -----

  ncks -O \
    -v q_reach \
    -d seg,"$i_summa","$i_summa" \
    -d method,0,0 \
    "$SUMMA_FILE" \
    "$summa_reach"

  # remove singleton method dimension
  ncwa -O \
    -a method \
    "$summa_reach" \
    "$summa_reach"

  # standalone output is float, so cast coupled output to float
  ncap2 -O \
    -s 'KWroutedRunoff=float(q_reach)' \
    "$summa_reach" \
    "$summa_reach"

  # compare the same-named variables
  ncdiff -O \
    -v KWroutedRunoff \
    "$mizu_reach" \
    "$summa_reach" \
    "$diff_reach"

  ncap2 -O \
    -s 'max_diff=max(abs(KWroutedRunoff))' \
    "$diff_reach" \
    "$summary_reach"

  reach_diff=$(ncks -H -C -s '%g' \
    -v max_diff \
    "$summary_reach")

  printf "  reachID=%-6s standalone_index=%-3s coupled_index=%-3s max abs diff=%s\n" \
       "$reach_id" "$i_mizu" "$i_summa" "$reach_diff"

  if awk "BEGIN {exit !($reach_diff > $max_diff)}"; then
    max_diff=$reach_diff
  fi

  if ! awk "BEGIN {exit !($reach_diff == 0)}"; then
    flow_result="FAIL"
  fi

  n_reaches=$((n_reaches + 1))

  if awk "BEGIN {exit !($reach_diff == 0)}"; then
    n_pass=$((n_pass + 1))
  else
    n_fail=$((n_fail + 1))
    flow_result="FAIL"
  fi
  
  if awk "BEGIN {exit !($reach_diff > $max_diff)}"; then
    max_diff=$reach_diff
    max_diff_reach=$reach_id
  fi

done

# ==============================================================================
# Regression summary
# ==============================================================================

echo
echo "======================================================================"
echo "mizuRoute coupling regression summary"
echo "======================================================================"

printf "%-32s %12s\n" "Reaches compared:"        "$n_reaches"
printf "%-32s %12s\n" "Exact matches:"           "$n_pass"
printf "%-32s %12s\n" "Reaches with differences:" "$n_fail"
printf "%-32s %12s\n" "Maximum absolute diff:"   "$max_diff"

if [[ -n "$max_diff_reach" ]]; then
  printf "%-32s %12s\n" "Reach with max diff:" "$max_diff_reach"
fi

echo "----------------------------------------------------------------------"

if [[ "$flow_result" == "PASS" ]]; then
  echo "PASS: all routed streamflow values match exactly"
else
  echo "FAIL: routed streamflow differences detected"
fi

echo "======================================================================"
echo


