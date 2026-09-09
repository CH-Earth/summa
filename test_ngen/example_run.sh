#!/bin/bash
# Run the gauge_01073000 SUMMA-BMI + t-route example.
# Run this from the main ngen directory:  ./extern/summa/summa/test_ngen/example_run.sh

# --- Python environment -------------------------------------------------------
# t-route (nwm_routing) must be importable, and ngen's embedded interpreter must
# see only ONE set of site-packages.  Activate the ngen env and pin VIRTUAL_ENV
# to it: ngen globs "**/site-packages/" under $VIRTUAL_ENV, so a stale value
# pointing at the conda root (as VS Code's Python extension sometimes sets) drags
# every other env's packages onto sys.path and crashes the interpreter.
: "${NGEN_CONDA_ENV:=ngen}"
if [ -f "${HOME}/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/opt/anaconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  . "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)" || true
fi
if command -v conda >/dev/null 2>&1; then
  conda activate "${NGEN_CONDA_ENV}" || true
fi
: "${NGEN_PYTHON:=$(command -v python || echo /usr/bin/python3)}"
export VIRTUAL_ENV="$("${NGEN_PYTHON}" -c 'import sys; print(sys.prefix)' 2>/dev/null || echo "${CONDA_PREFIX}")"
"${NGEN_PYTHON}" -c 'import nwm_routing' 2>/dev/null || {
  echo "ERROR: 'nwm_routing' not importable with ${NGEN_PYTHON}."
  echo "       Build t-route into the '${NGEN_CONDA_ENV}' env: cd extern/t-route && ./compiler_mac.sh"
  exit 1
}
# ---------------------------------------------------------------------------

#./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/gauge_01073000/settings/example_realization_config_w_summa_bmi_routing.json

./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./extern/summa/summa/test_ngen/gauge_01073000/settings/example_realization_config_w_summa_bmi.json
"${NGEN_PYTHON}" -m nwm_routing -V4 -f  ./test/data/routing/ngen_routing_config_unit_test.yaml

#./cmake_build/ngen ./test/data/routing/gauge_01073000.gpkg '' ./test/data/routing/gauge_01073000.gpkg '' ./data/gauge_01073000/example_bmi_multi_realization_config_w_routing.json

rm -f ./test/data/routing/*.parquet
