#!/bin/bash
#
# Create a mizuRoute runoff-remapping file for a lumped SUMMA simulation.
#
# The SUMMA forcing file must contain a single HRU. Its hruId is mapped to
# every routing HRU in the mizuRoute topology file with a weight of 1.0.
#
#                       SUMMA HRU
#                           |
#              +------------+------------+
#              |            |            |
#              v            v            v
#         routing HRU 1  routing HRU 2  routing HRU N
#
# Example
#
# forcing=/path/to/summa_forcing.nc
# topology=/path/to/topology.nc
# remapping=/path/to/lumped_to_hru.nc
#
# utils/remapping/create_lumped_to_hru_mapping.sh \
#     "$forcing" \
#     "$topology" \
#     "$remapping"
#
# Requires NCO:
#   ncks
#   ncap2
#   ncrename
#

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 forcing.nc topology.nc output.nc"
    exit 1
fi

forcing=$1
topology=$2
output=$3

# ----- check input files -----

[[ -f "$forcing" ]] || {
    echo "ERROR: forcing file not found: $forcing"
    exit 1
}

[[ -f "$topology" ]] || {
    echo "ERROR: topology file not found: $topology"
    exit 1
}

# ----- get the single SUMMA HRU ID -----

summa_hru_id=$(ncks -H -C -s '%d\n' -v hruId "$forcing" || true)

if [[ $(printf "%s\n" "$summa_hru_id" | wc -l | tr -d ' ') -ne 1 ]]; then
    echo "ERROR: forcing file must contain exactly one hruId"
    exit 1
fi

echo
echo "Creating lumped-to-HRU runoff mapping"
echo "  SUMMA HRU ID : $summa_hru_id"
echo "  topology     : $topology"
echo "  output       : $output"
echo

# ----- temporary file -----

output_dir=$(dirname "$output")
output_name=$(basename "$output")
work_dir="${output_dir}/work"
work_file="${work_dir}/${output_name%.nc}_work.nc"

mkdir -p "$work_dir"

rm -f "$work_file"
rm -f "$output"

# ----- extract routing HRU IDs -----

ncks -O -h -C \
    -v hruId \
    "$topology" \
    "$work_file"

# ----- rename the dimension -----
#
# At this point:
#     hruId(hru)
#
# becomes:
#     hruId(polyid)

ncrename -O -d hru,polyid "$work_file"

# ----- create remapping variables -----

ncap2 -O -h -s "
    defdim(\"data\",\$polyid.size);

    polyid[polyid]    = int(hruId);
    nOverlaps[polyid] = 1;

    weight[data]      = 1.0;
    qhru_id[data]     = ${summa_hru_id};
" \
"$work_file" \
"$output"

# ----- retain only remapping variables -----

ncks -O -h \
    -v polyid,nOverlaps,weight,qhru_id \
    "$output" \
    "$output"

# ----- clean up -----

rm -f "$work_file"

echo "Created: $output"
