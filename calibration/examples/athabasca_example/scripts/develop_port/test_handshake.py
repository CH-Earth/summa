"""Step 1 test: SPOTPY -> Python (ctypes) -> Fortran C wrapper, with no SUMMA run.

Part A  sends three known parameter vectors to summa_echo and checks that Fortran returns
        checksum = sum(i * value_i) exactly (order and double precision intact).
Part B  runs a short SPOTPY DDS (200 samples) whose 'model' is summa_echo: SPOTPY samples the
        15 Athabasca parameters within their bounds, each vector goes through the wrapper, and
        the objective is how close the Fortran checksum is to a target. SPOTPY improving the
        objective shows values flow both ways through the full SPOTPY -> C path.
Usage:  SUMMA_LIB=<path to libsumma.so> python test_handshake.py
"""
import csv, ctypes, os, sys
from pathlib import Path
import numpy as np

lib_path = os.environ.get("SUMMA_LIB")
assert lib_path and Path(lib_path).exists(), "set SUMMA_LIB to the develop libsumma.so"
lib = ctypes.CDLL(lib_path)
lib.summa_echo.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.c_int]
lib.summa_echo.restype = None

params = list(csv.DictReader(open(Path(__file__).resolve().parents[1] / "athabasca_params.csv")))
names = [p["name"] for p in params]
low = np.array([float(p["low"]) for p in params]); high = np.array([float(p["high"]) for p in params])
w = np.arange(1, len(names) + 1, dtype=float)

def echo(x, verbose=0):
    arr = (ctypes.c_double * len(x))(*x); out = ctypes.c_double(); err = ctypes.c_int()
    lib.summa_echo(arr, len(x), ctypes.byref(out), ctypes.byref(err), verbose)
    sys.stdout.flush()
    return out.value, err.value

print("== Part A: known vectors ==")
ok = True
for label, x in (("defaults", [float(p["default"]) for p in params]), ("lower bounds", low), ("upper bounds", high)):
    got, err = echo(list(x), verbose=1 if label == "defaults" else 0); want = float(np.dot(w, x))
    good = err == 0 and abs(got - want) <= 1e-12 * max(1.0, abs(want))
    ok &= good
    print(f"{label:12s} fortran {got:.12e}  python {want:.12e}  {'PASS' if good else 'FAIL'}")
_, err = echo([1.0] * 14)
print(f"wrong length (14) -> err={err}  {'PASS' if err == 1 else 'FAIL'}"); ok &= (err == 1)

print("\n== Part B: SPOTPY DDS through the wrapper ==")
import spotpy
target = float(np.dot(w, (low + high) / 2))

class Setup:
    def __init__(self):
        self.params = [spotpy.parameter.Uniform(n, low=float(a), high=float(b)) for n, a, b in zip(names, low, high)]
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    def simulation(self, vector):
        val, err = echo([float(v) for v in vector])
        return [val if err == 0 else np.nan]
    def evaluation(self):
        return [target]
    def objectivefunction(self, simulation, evaluation):
        return -abs(simulation[0] - evaluation[0]) / abs(evaluation[0])   # 0 is perfect

import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):                          # hide SPOTPY's progress chatter
    sampler = spotpy.algorithms.dds(Setup(), dbname="handshake_dds", dbformat="ram", save_sim=False)
    sampler.sample(200)
res = sampler.getdata()
like = np.array(res["like1"])
print(f"samples {len(like)} | first-10 best {like[:10].max():.4f} | final best {like.max():.6f} (0 = target hit)")
improved = like.max() > like[:10].max()
print(f"SPOTPY improved through the wrapper: {'PASS' if improved else 'FAIL'}")
print("\nHANDSHAKE", "PASS" if (ok and improved) else "FAIL")
