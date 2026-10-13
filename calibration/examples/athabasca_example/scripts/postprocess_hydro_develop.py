"""Post-process full-period develop runs (run_best_hydrograph_develop.py) without re-running the model.
Reads coupled mizuRoute streamflow Q_reach(time, seg, method) from SUMMA's *_timestep.nc at segment 222,
makes daily means (24 h ending at 00:00, labelled at the end, as J() and the period-end obs),
and writes the daily CSV, the KGE summary and the hydrograph for each point.
Usage: python postprocess_hydro_develop.py <tag> [point ...]
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
from netCDF4 import Dataset, num2date
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

basin = Path(__file__).resolve().parent.parent
tag = sys.argv[1]
points = sys.argv[2:] or ["develop_best", "spotpy_best", "symfluence_best"]
SEG = 222
CAL, EVAL = ("2016-01-01", "2019-12-31"), ("2020-01-01", "2024-12-31")

def times(ds):
    t = ds["time"]
    return pd.to_datetime([str(x) for x in num2date(t[:], t.units, getattr(t, "calendar", "standard"))], format="ISO8601").round("min")

d = Dataset(basin / "observations" / "Athabasca_daily_flow_observations_periodend.nc")
obs = pd.Series(np.ma.filled(np.ma.asarray(d["q_obs"][:], dtype=float), np.nan), index=times(d)); d.close()
obs.index = obs.index.normalize()

def kge(out, period):
    x = out.loc[period[0]:period[1]].dropna()
    r = np.corrcoef(x.q_sim, x.q_obs)[0, 1]
    a, b = x.q_sim.std() / x.q_obs.std(), x.q_sim.mean() / x.q_obs.mean()
    return 1 - np.sqrt((r - 1)**2 + (a - 1)**2 + (b - 1)**2), r, a, b, len(x)

for point in points:
    work = basin / "output" / f"hydro_develop_{tag}" / point
    ts = sorted(work.glob("*_timestep.nc"))
    assert ts, f"no *_timestep.nc in {work}"
    ds = Dataset(ts[0])
    ids = np.asarray(ds["seg"][:]).ravel()
    assert SEG in ids, f"segment {SEG} not in {ts[0].name}"
    j = int(np.where(ids == SEG)[0][0])
    q = pd.Series(np.ma.filled(np.ma.asarray(ds["Q_reach"][:, j, 0], dtype=float), np.nan), index=times(ds))
    ds.close()
    qd = q[~q.index.duplicated()].sort_index().resample("D", closed="right", label="right").mean()
    out = pd.DataFrame({"q_sim": qd, "q_obs": obs.reindex(qd.index)})
    out.index.name = "date"
    out.to_csv(work / f"seg222_daily_{point}.csv", float_format="%.3f")

    kj = None   # J() value printed by the run itself
    for f in basin.glob(f"athabasca_hydro_dev_{tag}_*.out"):
        m = re.search(rf"RESULT {point}\s+J\(\) KGE \(calibration 2016-2019\) = ([0-9.]+)", f.read_text(errors="ignore"))
        if m: kj = float(m.group(1))
    k = {name: kge(out, per) for name, per in (("calibration", CAL), ("evaluation", EVAL))}
    lines = [f"point {point}   J() calibration KGE {kj}"]
    for name, per in (("calibration", CAL), ("evaluation", EVAL)):
        v, r, a, b, n = k[name]
        lines.append(f"{name:12s} {per[0]} to {per[1]}  KGE {v:.4f}  r {r:.3f}  alpha {a:.3f}  beta {b:.3f}  days {n}")
    (work / f"kge_{point}.txt").write_text("\n".join(lines) + "\n")
    for l in lines: print("RESULT " + l)

    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.axvspan(pd.Timestamp(CAL[0]), pd.Timestamp(CAL[1]), color="tab:green", alpha=0.08, label=f"Calibration 2016-2019, KGE {k['calibration'][0]:.3f}")
    ax.axvspan(pd.Timestamp(EVAL[0]), pd.Timestamp(EVAL[1]), color="tab:orange", alpha=0.08, label=f"Evaluation 2020-2024, KGE {k['evaluation'][0]:.3f}")
    ax.plot(out.index, out.q_obs, color="black", lw=0.8, label="Observed")
    ax.plot(out.index, out.q_sim, color="tab:blue", lw=0.8, label=f"Simulated ({point})")
    ax.set_ylabel("Discharge (m³/s)")
    ax.set_title(f"Athabasca, outlet segment {SEG}: {point} (develop build, coupled mizuRoute KW)")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(work / f"hydrograph_{point}.png", dpi=150); plt.close(fig)
    print(f"  wrote {work.name}/seg222_daily_{point}.csv, kge_{point}.txt, hydrograph_{point}.png")
