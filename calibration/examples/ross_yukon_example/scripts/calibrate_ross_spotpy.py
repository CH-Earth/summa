"""SPOTPY calibration of Ross with J(params).
Usage: python calibrate_ross_spotpy.py <dds|sceua> <repetitions> [tag]"""
import csv, os, sys, time
from pathlib import Path
import numpy as np
import spotpy

ross = Path(__file__).resolve().parent.parent
os.environ["SUMMA_MASTER"] = str(ross / "settings" / "SUMMA" / "fileManager.txt")
os.environ["SUMMA_CONFIG"] = str(ross / "settings" / "ross_config.toml")
sys.path.insert(0, str(ross.parent / "calibration"))
from summa_python_interface import J, PARAM_NAMES

ALGO = sys.argv[1] if len(sys.argv) > 1 else "dds"
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
TAG = sys.argv[3] if len(sys.argv) > 3 else time.strftime("%Y%m%d_%H%M")
MAXIMIZE = {"dds": True, "sceua": False}[ALGO]   # verified by the direction test
PENALTY = -1e6                                    # same value J() returns for a crashed run

with open(ross / "scripts" / "ross_params.csv") as f:
    rows = {r["name"]: r for r in csv.DictReader(f)}
assert list(rows) == PARAM_NAMES, "ross_params.csv order must match PARAM_NAMES"

import json, shutil
np.random.seed(42)                                # fixed seed, as used by the group
out = ross / "results" / f"{ALGO}_{TAG}"
out.mkdir(parents=True, exist_ok=True)
shutil.copy(ross / "settings" / "ross_config.toml", out)
shutil.copy(ross / "scripts" / "ross_params.csv", out)
log = open(out / "evaluations.csv", "w", buffering=1)
log.write("n,kge,seconds," + ",".join(PARAM_NAMES) + "\n")

class RossSetup:
    def __init__(self):
        self.params = [spotpy.parameter.Uniform(n, low=float(r["low"]), high=float(r["high"]),
                                                optguess=float(r["default"])) for n, r in rows.items()]
        self.n, self.best, self.best_x = 0, -np.inf, None
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    def simulation(self, x):
        x, t0 = [float(v) for v in x], time.time()
        p = dict(zip(PARAM_NAMES, x))
        kge = PENALTY if p["heightCanopyTop"] <= p["heightCanopyBottom"] else J(x)
        self.n += 1
        if kge > self.best:
            self.best, self.best_x = kge, x
        log.write(f"{self.n},{kge},{time.time() - t0:.1f}," + ",".join(map(str, x)) + "\n")
        print(f"PROGRESS {ALGO} eval {self.n}: KGE={kge:.4f}  best={self.best:.4f}", flush=True)
        return [kge]
    def evaluation(self):
        return [0.0]
    def objectivefunction(self, simulation, evaluation, params=None):
        return simulation[0] if MAXIMIZE else -simulation[0]

setup = RossSetup()
db = str(out / "spotpy_db")
if ALGO == "sceua":
    sampler = spotpy.algorithms.sceua(setup, dbname=db, dbformat="csv")
    sampler.sample(REPS, ngs=4)                   # 4 complexes, as NUMBER_OF_COMPLEXES in config_ross.yaml
else:
    sampler = spotpy.algorithms.dds(setup, dbname=db, dbformat="csv")
    x0 = np.array([float(r["default"]) for r in rows.values()])
    sampler.sample(REPS, trials=1, x_initial=x0)  # start from the defaults, like SYMFLUENCE
json.dump({"algorithm": ALGO, "evaluations": setup.n, "best_kge": setup.best,
           "best_params": dict(zip(PARAM_NAMES, setup.best_x or []))},
          open(out / "best_params.json", "w"), indent=2)
print(f"DONE {ALGO}: {setup.n} evaluations, best KGE {setup.best:.4f}  -> {out}")
