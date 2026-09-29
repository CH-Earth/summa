"""One evaluation of J() on Ross, with SYMFLUENCE's best parameters."""
import os, sys, time
from pathlib import Path

ross = Path(__file__).resolve().parent.parent
os.environ["SUMMA_MASTER"] = str(ross / "settings" / "SUMMA" / "fileManager.txt")
os.environ["SUMMA_CONFIG"] = str(ross / "settings" / "ross_config.toml")
sys.path.insert(0, str(ross.parent / "calibration"))
from summa_python_interface import J, PARAM_NAMES

best = {"k_soil": 0.00941684, "theta_sat": 0.613238, "aquiferBaseflowExp": 1.52518,
        "aquiferBaseflowRate": 3.05587e-07, "qSurfScale": 3.34253, "summerLAI": 9.31163,
        "frozenPrecipMultip": 0.987852, "Fcapil": 0.0598626, "tempCritRain": 273.647,
        "heightCanopyTop": 18.2841, "heightCanopyBottom": 4.93957, "windReductionParam": 0.912355,
        "vGn_n": 3.55503, "routingGammaScale": 23435.2, "routingGammaShape": 3.6619}

t0 = time.time()
kge = J([best[n] for n in PARAM_NAMES])
print(f"KGE = {kge:.4f}   ({time.time() - t0:.1f} s)")
