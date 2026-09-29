"""
summa_python_interface.py

Thin ctypes wrapper around libsumma.dylib (built from summa_c_api.f90).
Exposes a plain J(params) -> float function that any Python optimizer
(scipy.optimize, SPOTPY, etc.) can call directly.

Usage:
    from summa_python_interface import J, PARAM_NAMES
    value = J([0.5, 0.3, 2.1])   # order must match PARAM_NAMES
"""

import ctypes
from pathlib import Path

# Must match PARAM_NAMES in summa_c_api.f90, in the same order
PARAM_NAMES = ["k_soil", "theta_sat", "aquiferBaseflowExp",
               "aquiferBaseflowRate", "qSurfScale", "summerLAI", "frozenPrecipMultip",
               "Fcapil", "tempCritRain", "heightCanopyTop", "heightCanopyBottom",
               "windReductionParam", "vGn_n", "routingGammaScale", "routingGammaShape"]
N_PARAMS = len(PARAM_NAMES)

# Adjust extension per platform: .dylib (mac), .so (linux), .dll (windows)
LIB_PATH = Path(__file__).parent / "libsumma.dylib"

# Repo root, so master/config paths work regardless of cwd
REPO_ROOT = Path(__file__).parent.parent
import os
MASTER_FILE = Path(os.environ.get("SUMMA_MASTER", REPO_ROOT / "test_coupled" / "settings" / "summa_fileManager.txt"))
CONFIG_FILE = Path(os.environ.get("SUMMA_CONFIG", REPO_ROOT / "test_coupled" / "settings" / "summa_config_test.toml"))

_lib = ctypes.CDLL(str(LIB_PATH))

# --- gfortran argv fix ---
# gfortran's command-line intrinsics read argv via _gfortran_set_args,
# which normally gets populated from a Fortran PROGRAM's real argc/argv.
# Since Python (not a Fortran/C main) is the process entry point here,
# we must call this manually with a synthetic argv before summa_evaluate().
_lib._gfortran_set_args.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
_lib._gfortran_set_args.restype = None

_argv_list = [
    b"summa_evaluate",
    b"-m", str(MASTER_FILE).encode(),
    b"-c", str(CONFIG_FILE).encode(),
    b"-s", b"_pytest1",
]
_argc = len(_argv_list)
_argv_arr = (ctypes.c_char_p * _argc)(*_argv_list)
_lib._gfortran_set_args(_argc, _argv_arr)
# --- end argv fix ---

_lib.summa_evaluate.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # param_values
    ctypes.c_int,                     # n
    ctypes.POINTER(ctypes.c_double),  # objective (out)
    ctypes.POINTER(ctypes.c_int),     # err (out)
]
_lib.summa_evaluate.restype = None


def J(params, penalty=-1e6):
    """
    Evaluate the SUMMA objective function for a given parameter set.

    params  : sequence of floats, length N_PARAMS, in PARAM_NAMES order
    penalty : value returned if SUMMA reports an error (e.g. non-physical
              parameter combination causing a model crash) — keeps
              optimizers from stalling on failed runs. Tune sign/magnitude
              once you know whether the objective is minimized or maximized.

    Returns: float objective value
    """
    if len(params) != N_PARAMS:
        raise ValueError(f"Expected {N_PARAMS} parameters, got {len(params)}")

    values_arr = (ctypes.c_double * N_PARAMS)(*params)
    objective  = ctypes.c_double()
    err        = ctypes.c_int()

    _lib.summa_evaluate(values_arr, N_PARAMS,
                        ctypes.byref(objective), ctypes.byref(err))

    if err.value != 0:
        print(f"[WARN] SUMMA evaluation failed (err={err.value}) "
              f"for params={list(params)}; returning penalty value.")
        return penalty

    return objective.value


if __name__ == "__main__":
    # Same params as today's validated C-call baseline:
    # k_soil=7.5e-6, theta_sat=0.55, vGn_n=1.3 -> objective = -10.2181321285
    test_params = [7.5e-06, 0.55, 1.3]
    print(f"J({test_params}) = {J(test_params)}")