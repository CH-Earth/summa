#!/usr/bin/env bash

set -euo pipefail

ROOT="$HOME/data/century/test/exp01/domain"

for case_dir in "${ROOT}"/*; do

  [[ -d "${case_dir}" ]] || continue

  case_name=$(basename "${case_dir}")
  forcing_file="${case_dir}/summa_forcing/${case_name}_em_earth_distributed.nc"

  if [[ ! -f "${forcing_file}" ]]; then
    echo "${case_name}: forcing file not found"
    continue
  fi

  first_time=$(ncks -H --cal -v time "${forcing_file}" | \
               sed -n 's/.*time = "\([^"]*\)".*/\1/p' | head -1)

  last_time=$(ncks -H --cal -v time "${forcing_file}" | \
              grep -o '"[^"]*"' | tail -1 | tr -d '"')

  printf '%-15s  %s  to  %s\n' \
    "${case_name}" "${first_time}" "${last_time}"

done
