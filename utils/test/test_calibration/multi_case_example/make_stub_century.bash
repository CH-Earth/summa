#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------
# Build a stub Century dataset, so the multi-case calibration path can be run before the
# real dataset is available.
#
# The real Century dataset is ~115 basins that are not in this repository. This assembles
# the same *directory layout* that summa_config_template.toml expects, with a few cases
# whose inputs all point at the one domain that is in the repository, the Bow at Banff
# (utils/test/test_mizuroute/bow_real_data). Every case is therefore the same basin under
# a different name -- that is fine for exercising the manifest, the case distribution and
# the per-case output, which is what has no other test.
#
# Files are symlinked, not copied, so the stub costs almost nothing on disk.
#
# It also writes a stub manifest and template beside the data, derived from the committed
# manifest_century.toml and summa_config_template.toml but retargeted at the stub root and
# shortened to a period the bundled forcing covers. The committed originals are left alone
# as the reference for the real dataset.
#
# Usage:  ./make_stub_century.bash [root] [n_cases] [cases_per_node]
#           root            where to build it  (default: a stub_century directory beside this script)
#           n_cases         how many cases     (default: 3)
#           cases_per_node  cases run at once  (default: 1)
#
# NOTE: the driver requires the MPI rank count to be divisible by cases_per_node, and each case
#       group needs at least two ranks (one coordinates, the rest evaluate samples). So run with
#       -np equal to 2 * cases_per_node. The script prints the right command when it finishes.
# ---------------------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMMA_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BOW="${SUMMA_ROOT}/utils/test/test_mizuroute/bow_real_data"

ROOT=${1:-"${SCRIPT_DIR}/stub_century"}
N_CASES=${2:-3}
CASES_PER_NODE=${3:-1}

if [ "${CASES_PER_NODE}" -gt "${N_CASES}" ]; then
  echo "ERROR: cases_per_node (${CASES_PER_NODE}) cannot exceed n_cases (${N_CASES})"
  exit 1
fi

if [ ! -d "${BOW}" ]; then
  echo "ERROR: bundled Bow domain not found at ${BOW}"
  exit 1
fi

# the real dataset's first cases, so the stub uses names that will still mean something
ALL_CASES=(CAN_05BB001 CAN_01AD003 CAN_01AM001 CAN_01BJ010 CAN_01DG003 CAN_01FB001)
CASES=("${ALL_CASES[@]:0:${N_CASES}}")

FORCING="CAN_05BB001_em_earth_distributed_1980-1990.nc"

echo "Building stub Century dataset"
echo "  root:    ${ROOT}"
echo "  cases:   ${CASES[*]}"
echo "  source:  ${BOW}"
echo

rm -rf "${ROOT}"

# ----- common configuration files, shared by every case ---------------------------------
COMMON="${ROOT}/models/summa/test_coupled/case/settings"
mkdir -p "${COMMON}"
for f in modelDecisions.txt outputControl_short.txt trialParams.nc \
         localParamInfo.txt basinParamInfo.txt \
         TBL_VEGPARM.TBL TBL_SOILPARM.TBL TBL_GENPARM.TBL TBL_MPTABLE.TBL; do
  ln -s "${BOW}/summa_inputs/${f}" "${COMMON}/${f}"
done

# ----- per-case inputs ------------------------------------------------------------------
OBS_DIR="${ROOT}/data/camels-spat/observations/obs-daily"
mkdir -p "${OBS_DIR}"

for case_name in "${CASES[@]}"; do
  case_dir="${ROOT}/data/century/test/exp01/domain/${case_name}"
  mkdir -p "${case_dir}"/{summa_inputs,summa_forcing,summa_state,mizuroute_inputs,work}

  ln -s "${BOW}/summa_inputs/attributes.nc" "${case_dir}/summa_inputs/attributes.nc"
  ln -s "${BOW}/summa_inputs/${FORCING}"    "${case_dir}/summa_forcing/${FORCING}"
  ln -s "${BOW}/summa_inputs/coldState.nc"  "${case_dir}/summa_state/coldState.nc"

  # the forcing list is read relative to settings_path, so it names the file only
  echo "${FORCING}" > "${case_dir}/summa_inputs/forcingFileList.txt"

  for f in topology.nc lumped_to_hru.nc mizuroute.param; do
    ln -s "${BOW}/mizuroute_inputs/${f}" "${case_dir}/mizuroute_inputs/${f}"
  done

  # observations are looked up by case name
  ln -s "${BOW}/mizuroute_inputs/CAN_05BB001_daily_flow_observations.nc" \
        "${OBS_DIR}/${case_name}_daily_flow_observations.nc"

  echo "  created ${case_name}"
done

# ----- stub manifest and template --------------------------------------------------------
case_list=$(printf '"%s", ' "${CASES[@]}"); case_list="${case_list%, }"

sed -e "s|^template_path = .*|template_path = \"${ROOT}/\"|" \
    -e "s|^template_file = .*|template_file = \"stub_template.toml\"|" \
    -e "s|^cases_per_node = .*|cases_per_node = ${CASES_PER_NODE}|" \
    "${SCRIPT_DIR}/manifest_century.toml" \
  | awk -v list="${case_list}" '
      /^case_names = \[/ { print "case_names = [ " list " ]"; skip=1; next }
      skip && /\]/       { skip=0; next }
      skip               { next }
      { print }
    ' > "${ROOT}/stub_manifest.toml"

# retarget the template at the stub root, and shorten to what the bundled forcing covers
# (1980-01-01 to 1990-01-01), leaving 1980 as the spinup year
sed -e "s|^home_path = .*|home_path = \"${ROOT}\"|" \
    -e "s|^start_time = .*|start_time = \"1981-01-01 00:00\"|" \
    -e "s|^end_time   = .*|end_time   = \"1981-12-31 23:00\"|" \
    -e "s|^start_date     = .*|start_date     = \"1981-01-01\"|" \
    -e "s|^end_date       = .*|end_date       = \"1981-12-31\"|" \
    -e "s|^write_aligned  = .*|write_aligned  = false\nn_samples      = 3|" \
    "${SCRIPT_DIR}/summa_config_template.toml" > "${ROOT}/stub_template.toml"

echo
echo "wrote ${ROOT}/stub_manifest.toml"
echo "wrote ${ROOT}/stub_template.toml"
echo
NP=$((CASES_PER_NODE * 2))
echo "run it with (-np must be divisible by cases_per_node = ${CASES_PER_NODE}):"
echo "  mpirun -np ${NP} ${SUMMA_ROOT}/bin/summa_sundials_mizuroute_opt.exe \\"
echo "         --manifest ${ROOT}/stub_manifest.toml"
