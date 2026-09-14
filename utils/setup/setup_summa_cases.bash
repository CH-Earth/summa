#!/usr/bin/env bash

# Prepare SUMMA multi-case input directories.
#
# For each case listed in CASE_LIST, this script copies the required SUMMA
# and mizuRoute inputs from the existing Century dataset into the directory
# structure used for multi-case simulations. It also creates the mapping
# between the SUMMA HRU and the mizuRoute HRUs.

set -e


# ============================================================================
# User-configurable settings
# ============================================================================

# Root directory containing the existing Century model cases.
OLD_ROOT="$HOME/data/century/MM"

# Root directory for the new multi-case experiment.
NEW_ROOT="$HOME/data/century/test/exp01"

# File containing the list of cases to prepare, one case name per line.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_LIST="${SCRIPT_DIR}/century_cases.txt"


# ============================================================================
# Script configuration
# ============================================================================

# Locate the SUMMA utilities directory relative to this script.
UTILS_DIR="$(dirname "${SCRIPT_DIR}")"


# ============================================================================
# Prepare cases
# ============================================================================

while IFS= read -r case_name; do

  # Skip empty and comment lines.
  [[ -z "${case_name}" ]] && continue
  [[ "${case_name}" == \#* ]] && continue

  # Define the source and destination directories for this case.
  old_case="${OLD_ROOT}/${case_name}/7-FA_mod_IC_newSUMMA_inf_GA_bsflwParams_MP/common_inputs"
  new_case="${NEW_ROOT}/domain/${case_name}"

  echo "Creating case: ${case_name}"

  # Create the directory structure for the new case.
  mkdir -p \
    "${new_case}/summa_state" \
    "${new_case}/summa_inputs" \
    "${new_case}/summa_forcing" \
    "${new_case}/mizuroute_inputs" \
    "${new_case}/work"

  # Copy SUMMA initial conditions.
  cp "${old_case}/summa/coldState.nc" \
     "${new_case}/summa_state/"

  # Copy SUMMA model inputs.
  cp "${old_case}/summa/attributes.nc" \
     "${new_case}/summa_inputs/"
  cp "${old_case}/summa/forcingFileList.txt" \
     "${new_case}/summa_inputs/"
  cp "${old_case}/summa/trialParams.priori.nc" \
     "${new_case}/summa_inputs/"

  # Copy SUMMA meteorological forcing.
  cp "${old_case}/summa/${case_name}_em_earth_distributed.nc" \
     "${new_case}/summa_forcing/"

  # Copy mizuRoute inputs.
  cp "${old_case}/mizuroute/topology.nc" \
     "${new_case}/mizuroute_inputs/"
  cp "${old_case}/mizuroute/mizuroute.param" \
     "${new_case}/mizuroute_inputs/"

  # Create the mapping between the SUMMA HRU and mizuRoute HRUs.
  "${UTILS_DIR}/remapping/create_lumped_to_hru_mapping.sh" -q \
    "${new_case}/summa_forcing/${case_name}_em_earth_distributed.nc" \
    "${new_case}/mizuroute_inputs/topology.nc" \
    "${new_case}/mizuroute_inputs/lumped_to_hru.nc"

done < "${CASE_LIST}"
