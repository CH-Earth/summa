"""Full-period run of one parameter set with the NCAR/summa develop build + coupled mizuRoute (KW),
full output (writeTS) and period-end observations.
Period: spin-up 2014-09-01, calibration 2016-2019, evaluation 2020-2024 (settings/SUMMA/fileManager_full.txt).
Usage: python run_best_hydrograph_develop.py <develop_best|spotpy_best|symfluence_best|default> [tag]
Output: output/hydro_develop_<tag>/<point>/
  seg222_daily_<point>.csv  date, q_sim, q_obs  (daily = 24 h ending at the date's 00:00, as J() and the obs file)
  kge_<point>.txt           KGE, r, alpha, beta for the calibration and evaluation periods
  hydrograph_<point>.png
"""
import csv, json, os, re, sys, time
from pathlib import Path

basin = Path(__file__).resolve().parent.parent
ROOT = basin.parent.parent
point = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "local"))
SEG = 222
CAL, EVAL = ("2016-01-01", "2019-12-31"), ("2020-01-01", "2024-12-31")

with open(basin / "scripts" / "athabasca_params.csv") as f:
    default = {r["name"]: float(r["default"]) for r in csv.DictReader(f)}
res = basin / "results"
develop_best = json.load(open(res / "asyncdds_mpi_develop_noTS_periodend_21217337" / "summary.json"))["best_params"]
spotpy_best = json.load(open(res / "asyncdds_mpi_21150407" / "summary.json"))["best_params"]
symfluence_best = {"k_soil": 0.0009428590809408767, "theta_sat": 0.39871409009340836,
    "aquiferBaseflowExp": 2.926994945210623, "aquiferBaseflowRate": 1.837279373512553e-06,
    "qSurfScale": 6.531566859453017, "summerLAI": 8.290895744495714,
    "frozenPrecipMultip": 0.5611296554864641, "Fcapil": 0.0053496416060497464,
    "tempCritRain": 272.16554891743357, "heightCanopyTop": 4.780397129075084,
    "heightCanopyBottom": 3.2272406576351074, "windReductionParam": 0.3599284468484365,
    "vGn_n": 1.2648314071806919, "routingGammaScale": 153.00534938182162,
    "routingGammaShape": 2.0046760039588944}
p = {"develop_best": develop_best, "spotpy_best": spotpy_best,
     "symfluence_best": symfluence_best, "default": default}[point]

work = basin / "output" / f"hydro_develop_{tag}" / point
work.mkdir(parents=True, exist_ok=True)
def fm_line(l):
    k = l.strip()
    if k.startswith("outputPath"):
        return f"outputPath '{work}/'"
    if k.startswith("outputControlFile"):
        return "outputControlFile 'outputControl_develop.txt'"   # develop's renamed output variable
    return l
fm = [fm_line(l) for l in (basin / "settings" / "SUMMA" / "fileManager_full.txt").read_text().splitlines()]
(work / "fileManager.txt").write_text("\n".join(fm) + "\n")
cfg = (basin / "settings" / "athabasca_config_develop_writeTS.toml").read_text()
cfg = re.sub(r'^hfabric_path\s*=.*$', f'hfabric_path = "{work}/"', cfg, flags=re.M)
cfg = re.sub(r'^write_Qreach\s*=.*$', 'write_Qreach  = true', cfg, flags=re.M)
(work / "config.toml").write_text(cfg)
if not (work / "topology.nc").exists():
    os.symlink(basin / "settings" / "mizuRoute" / "topology.nc", work / "topology.nc")
os.environ["SUMMA_MASTER"] = str(work / "fileManager.txt")
os.environ["SUMMA_CONFIG"] = str(work / "config.toml")
os.environ["SUMMA_LIB"] = str(ROOT / "summa-develop" / "build" / "cmake_build" / "libsumma.so")
sys.path.insert(0, str(basin / "scripts" / "develop_port"))
os.chdir(work)
from summa_python_interface import J, PARAM_NAMES

t0 = time.time()
kge_J = J([p[n] for n in PARAM_NAMES])
print(f"RESULT {point:16s} J() KGE (calibration 2016-2019) = {kge_J:.4f}   wall time {time.time()-t0:.0f} s", flush=True)

import numpy as np, pandas as pd
from netCDF4 import Dataset, num2date
def series(ds, values):
    t = ds["time"]
    idx = pd.to_datetime([str(x) for x in num2date(t[:], t.units, getattr(t, "calendar", "standard"))])
    return pd.Series(np.ma.filled(np.ma.asarray(values, dtype=float), np.nan), index=idx)

q = None
for nc in sorted(work.rglob("*.nc")):
    if nc.name == "topology.nc":
        continue
    try:
        ds = Dataset(nc)
    except Exception:
        continue
    names = [v for v in ds.variables if re.search(r"routedRunoff$", v, re.I)]
    idv = next((v for v in ("reachID", "segId", "seg_id") if v in ds.variables), None)
    if names and idv and "time" in ds.variables:
        ids = np.asarray(ds[idv][:]).ravel()
        if SEG in ids:
            var = ds[names[0]]
            axis = var.dimensions.index(ds[idv].dimensions[-1])
            s = series(ds, np.take(var[:], int(np.where(ids == SEG)[0][0]), axis=axis))
            q = s if q is None else pd.concat([q, s])
            print(f"read {names[0]} from {nc.name}", flush=True)
    ds.close()
assert q is not None, f"no routed-runoff output with segment {SEG} found under {work}"
q = q[~q.index.duplicated()].sort_index()
qd = q.resample("D", closed="right", label="right").mean()     # 24 h ending at 00:00, labelled at the end
d = Dataset(basin / "observations" / "Athabasca_daily_flow_observations_periodend.nc")
obs = series(d, d["q_obs"][:]); d.close()
obs.index = obs.index.normalize()
out = pd.DataFrame({"q_sim": qd, "q_obs": obs.reindex(qd.index)})
out.index.name = "date"
out.to_csv(work / f"seg222_daily_{point}.csv", float_format="%.3f")

def kge(period):
    x = out.loc[period[0]:period[1]].dropna()
    r = np.corrcoef(x.q_sim, x.q_obs)[0, 1]
    a, b = x.q_sim.std() / x.q_obs.std(), x.q_sim.mean() / x.q_obs.mean()
    return 1 - np.sqrt((r - 1)**2 + (a - 1)**2 + (b - 1)**2), r, a, b, len(x)
lines = [f"point {point}   J() calibration KGE {kge_J:.4f}"]
k = {}
for name, per in (("calibration", CAL), ("evaluation", EVAL)):
    k[name] = kge(per)
    v, r, a, b, n = k[name]
    lines.append(f"{name:12s} {per[0]} to {per[1]}  KGE {v:.4f}  r {r:.3f}  alpha {a:.3f}  beta {b:.3f}  days {n}")
(work / f"kge_{point}.txt").write_text("\n".join(lines) + "\n")
for l in lines: print("RESULT " + l, flush=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(13, 4.5))
ax.axvspan(pd.Timestamp(CAL[0]), pd.Timestamp(CAL[1]), color="tab:green", alpha=0.08, label=f"Calibration 2016-2019, KGE {k['calibration'][0]:.3f}")
ax.axvspan(pd.Timestamp(EVAL[0]), pd.Timestamp(EVAL[1]), color="tab:orange", alpha=0.08, label=f"Evaluation 2020-2024, KGE {k['evaluation'][0]:.3f}")
ax.plot(out.index, out.q_obs, color="black", lw=0.8, label="Observed")
ax.plot(out.index, out.q_sim, color="tab:blue", lw=0.8, label=f"Simulated ({point})")
ax.set_ylabel("Discharge (m³/s)")
ax.set_title(f"Athabasca, outlet segment {SEG}: {point} (develop build, coupled mizuRoute KW)")
ax.legend(loc="upper left", fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(work / f"hydrograph_{point}.png", dpi=150)
print(f"wrote {work}/seg222_daily_{point}.csv, kge_{point}.txt, hydrograph_{point}.png", flush=True)
