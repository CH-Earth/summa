"""Convert the Fortran DDS calibration file (<run>/output/*calibration.nc, written by
build/source/calibration/calibration_output.f90) into the same files the Python drivers write, so that
compare_athabasca_optimizers.py can read it:
  <run>/results/evaluations.csv          every sample: n, worker, kge, seconds, elapsed_s, best_so_far, parameters
  <run>/results/evaluations_by_time.csv  run, hours, kge, best_so_far_kge (in order of completion)
  <run>/results/summary.json             best KGE and parameters, wall time, runs per hour, ...
Times: hours are counted from --start (the Slurm job start, as printed in athabasca_fortran_dds_<jobid>.out),
so they include the one-year spinup; without --start, from the first dispatched sample (spinup excluded).
Usage (from summa-spotpy/my_work_athabasca):
  python scripts/fortran_dds/convert_fortran_dds_output.py fortran_dds/full --start "2026-10-24 22:03:25" """
import argparse, csv, json, math
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("run_dir")
ap.add_argument("--start", help='job start, local time as in the .out file, e.g. "2026-10-24 22:03:25"')
a = ap.parse_args()
run = Path(a.run_dir)
files = sorted((run / "output").glob("*calibration.nc"))
assert files, f"no *calibration.nc in {run / 'output'}"
f = files[0]

# ---- read (netCDF4, else xarray, else h5py) ----
def read_nc(path):
    try:
        import netCDF4
        with netCDF4.Dataset(path) as ds:
            ds.set_auto_mask(False)
            return {k: np.array(v[:]) for k, v in ds.variables.items()}, \
                   {k: {n: v.getncattr(n) for n in v.ncattrs()} for k, v in ds.variables.items()}
    except ImportError:
        pass
    try:
        import xarray as xr
        ds = xr.open_dataset(path, decode_times=False, mask_and_scale=False)
        return {k: ds[k].values for k in ds.variables}, {k: dict(ds[k].attrs) for k in ds.variables}
    except ImportError:
        pass
    import h5py
    with h5py.File(path, "r") as h:
        return {k: np.array(h[k]) for k in h.keys() if isinstance(h[k], h5py.Dataset)}, \
               {k: {n: (v.decode() if isinstance(v, bytes) else v) for n, v in h[k].attrs.items()}
                for k in h.keys() if isinstance(h[k], h5py.Dataset)}

var, att = read_nc(f)
params = [k for k, at in att.items() if "sampled" in at and str(at.get("sampled")).lower() in ("1", "true", "yes", ".true.")]
if not params:   # fall back: every per-sample double that is not bookkeeping
    params = [k for k in var if k not in ("sample", "worker_rank", "dispatch_time", "completion_time", "objective")
              and np.ndim(var[k]) == 1]
obj = np.asarray(var["objective"], dtype=float)
end = np.asarray(var["completion_time"]); beg = np.asarray(var["dispatch_time"])
if end.shape[0] == 8 and end.shape[-1] != 8:     # (time_component, sample) -> (sample, time_component)
    end, beg = end.T, beg.T

def to_dt(c):     # year month day UTC_offset(min) hour minute second millisecond, local time
    y, mo, d, _, h, mi, s, ms = (int(x) for x in c)
    if y <= 0 or y > 3000:
        return None
    return datetime(y, mo, d, h, mi, s) + timedelta(milliseconds=ms)

rows = []
for i in range(len(obj)):
    te, tb = to_dt(end[i]), to_dt(beg[i])
    if te is None or tb is None or not math.isfinite(obj[i]) or abs(obj[i]) > 1e30:
        continue                                   # sample never completed (run stopped early)
    rows.append({"n": int(var["sample"][i]) if "sample" in var else i + 1,
                 "worker": int(var["worker_rank"][i]) if "worker_rank" in var else -1,
                 "kge": float(obj[i]), "t_begin": tb, "t_end": te,
                 **{p: float(var[p][i]) for p in params}})
assert rows, "no completed samples found"
rows.sort(key=lambda r: r["t_end"])
t0 = datetime.strptime(a.start, "%Y-%m-%d %H:%M:%S") if a.start else min(r["t_begin"] for r in rows)

out = run / "results"; out.mkdir(exist_ok=True)
best = -math.inf; bestrow = None
with open(out / "evaluations.csv", "w", newline="") as fe, open(out / "evaluations_by_time.csv", "w") as ft:
    fe.write("n,worker,kge,seconds,elapsed_s,best_so_far," + ",".join(params) + "\n")
    ft.write("run,hours,kge,best_so_far_kge\n")
    for k, r in enumerate(rows, 1):
        if r["kge"] > -1e5 and r["kge"] > best:
            best, bestrow = r["kge"], r
        el = (r["t_end"] - t0).total_seconds()
        sec = (r["t_end"] - r["t_begin"]).total_seconds()
        fe.write(f"{r['n']},{r['worker']},{r['kge']},{sec:.1f},{el:.1f},{best}," + ",".join(repr(r[p]) for p in params) + "\n")
        ft.write(f"{k},{el / 3600:.4f},{r['kge']},{best}\n")
elapsed = (rows[-1]["t_end"] - t0).total_seconds()
secs = np.array([(r["t_end"] - r["t_begin"]).total_seconds() for r in rows])
summary = {"algorithm": "Fortran async DDS (summa_driver_opt)", "source_file": str(f),
           "evaluations": len(rows), "samples_in_file": int(len(obj)),
           "failed": int(sum(r["kge"] <= -1e5 for r in rows)),
           "elapsed_s": round(elapsed, 1), "clock_start": "job start (incl. spinup)" if a.start else "first dispatch (excl. spinup)",
           "runs_per_hour": round(len(rows) / elapsed * 3600, 1),
           "median_run_min": round(float(np.median(secs)) / 60, 2),
           "best_kge": best, "best_sample": bestrow["n"],
           "best_params": {p: bestrow[p] for p in params}}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps({k: v for k, v in summary.items() if k != "best_params"}, indent=2))
print(f"\nWrote {out}/evaluations.csv, evaluations_by_time.csv, summary.json")
