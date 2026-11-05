#!/bin/bash
# ---------------------------------------------------------------------------------------
# NSGA-II seeded from an earlier calibration: Wolverine Glacier.
#
# Runs test_calibration_wolverine.sh, then repeats it with half of the first generation seeded
# from that run's final population. Seeds must be its best members with the same objectives,
# and the rest of the generation must be new.
#
# Usage:  ./test_calibration_seed.sh [population] [generations] [n_ranks]
# ---------------------------------------------------------------------------------------
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK="${TEST_DIR}/wolverine_run"
POPULATION=${1:-4}
GENERATIONS=${2:-2}
N_RANKS=${3:-3}
CALIB=gulkana_wolverine_run0001_calibration.nc

"${TEST_DIR}/test_calibration_wolverine.sh" "${POPULATION}" "${GENERATIONS}" "${N_RANKS}"

SUMMA_EXE=""
for exe in summa_sundials_mizuroute_opt.exe summa_mizuroute_opt.exe; do
  if [ -x "${TEST_DIR}/../../../bin/${exe}" ]; then SUMMA_EXE="${TEST_DIR}/../../../bin/${exe}"; break; fi
done

echo
echo "Repeating with half of the first generation seeded..."
mv "${WORK}/output/${CALIB}" "${WORK}/earlier_calibration.nc"
sed -i.bak "s|^\[calibration.nsga2\]|[calibration.nsga2]\n\nseed_file = \"${WORK}/earlier_calibration.nc\"\nseed_fraction = 0.5|" "${WORK}/config.toml"
mpirun -np "${N_RANKS}" "${SUMMA_EXE}" -c "${WORK}/config.toml" > "${WORK}/run_seeded.log" 2>&1 || { tail -25 "${WORK}/run_seeded.log"; echo "FAIL: seeded calibration failed"; exit 1; }
grep -E "^NSGA-II" "${WORK}/run_seeded.log" | sed 's/^/  /' || true

python3 - "${WORK}/earlier_calibration.nc" "${WORK}/output/${CALIB}" "${POPULATION}" <<'PYEOF'
import sys
import numpy as np
import netCDF4 as nc

old, new, pop = nc.Dataset(sys.argv[1]), nc.Dataset(sys.argv[2]), int(sys.argv[3])
for d in (old, new):
    d.set_auto_mask(False)
nSeed = max(round(0.5 * pop), 1)
names = [v for v in old.variables if hasattr(old[v], "lower_bound") and old[v].sampled == 1]

# best first: lowest front, then largest crowding distance, without failed trials
member = old["population"][-1] - 1
rank = old["population_rank"][-1]
crowd = old["population_crowding"][-1]
failed = old["trial_failed"][:] == 1
order = [i for i in np.lexsort((-crowd, rank)) if not failed[member[i]]]
best = member[order][:nSeed]
if len(best) == 0:
    sys.exit("FAIL: the earlier run has no member that ran")

fail = False
def check(ok, what):
    global fail
    print(f"  {'ok  ' if ok else 'FAIL'} {what}")
    fail |= not ok

n = len(best)
seeded = np.array([new[v][:n] for v in names])
source = np.array([old[v][best] for v in names])
others = np.array([new[v][n:pop] for v in names])
check(np.array_equal(seeded, source), f"the first {n} of {pop} trials hold the earlier run's best members, bit for bit")
check(np.array_equal(new["objective"][:n], old["objective"][best]), "they score the same objectives")
earlier = np.array([old[v][:] for v in names])
check(not any(np.array_equal(c, r) for c in others.T for r in earlier.T), f"the other {pop - n} trials of the generation are new draws")
sys.exit(1 if fail else 0)
PYEOF

echo
echo "PASS"
