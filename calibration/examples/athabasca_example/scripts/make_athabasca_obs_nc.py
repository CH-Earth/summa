"""Write Athabasca observations in the same NetCDF layout as the Ross file J() already reads.
Uses my_work_ross/observations/Ross_daily_flow_observations.nc as the template (dimension, time
units, variable name and attributes), filled with the WSC daily discharge from SYMFLUENCE's CSV.
Usage: python make_athabasca_obs_nc.py
The series is regular daily (J() rejects irregular observation time steps)."""
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc

basin = Path(__file__).resolve().parent.parent
tmpl = basin.parent / "my_work_ross" / "observations" / "Ross_daily_flow_observations.nc"
src = basin / "observations" / "Athabasca_sedi_3_streamflow_processed.csv"
dst = basin / "observations" / "Athabasca_daily_flow_observations.nc"

df = pd.read_csv(src, usecols=[0, 1])
df.columns = ["date", "q"]
df["date"] = pd.to_datetime(df["date"], format="%Y/%m/%d")
df = df.dropna().sort_values("date")
# J() needs a regular daily time step: keep the model years only (forcing starts 2014-09-01) on a
# continuous daily axis, and fill the few missing days by linear interpolation (2020-05-01 to 05-04).
START, END = "2014-09-01", "2024-12-31"
s = df.set_index("date")["q"].loc[START:END]
full = pd.date_range(START, END, freq="D")
gaps = full.difference(s.index)
s = s.reindex(full).interpolate(limit_area="inside")
assert s.notna().all(), "missing values at the ends of the period"
df = pd.DataFrame({"date": full, "q": s.to_numpy()})
print(f"Interpolated {len(gaps)} missing days: {', '.join(gaps.strftime('%Y-%m-%d'))}")

with nc.Dataset(tmpl) as t, nc.Dataset(dst, "w") as o:
    q_t = t["q_obs"]
    assert q_t.ndim == 1, f"template q_obs has dims {q_t.dimensions}; adapt this script"
    tdim = q_t.dimensions[0]
    tvar = t[tdim] if tdim in t.variables else t["time"]
    o.setncatts({k: t.getncattr(k) for k in t.ncattrs()})
    o.title = "Athabasca_sedi_3 daily streamflow (WSC via SYMFLUENCE), 2014-09-01 to 2024-12-31, regular daily step, layout copied from Ross"
    o.gap_filling = "linear interpolation of " + ", ".join(gaps.strftime("%Y-%m-%d"))
    o.createDimension(tdim, None if t.dimensions[tdim].isunlimited() else len(df))
    tv = o.createVariable(tvar.name, tvar.dtype, (tdim,))
    tv.setncatts({k: tvar.getncattr(k) for k in tvar.ncattrs()})
    cal = getattr(tvar, "calendar", "standard")
    tv[:] = nc.date2num(df["date"].dt.to_pydatetime(), tvar.units, cal)
    fill = getattr(q_t, "_FillValue", None)
    qv = o.createVariable("q_obs", q_t.dtype, (tdim,), fill_value=fill)
    qv.setncatts({k: q_t.getncattr(k) for k in q_t.ncattrs() if k != "_FillValue"})
    qv[:] = df["q"].to_numpy(dtype=float)
    print(f"Template: {tmpl.name}  dims {q_t.dimensions}  time units '{tvar.units}'")
print(f"Wrote {dst}: {len(df)} days, {df['date'].iloc[0].date()} to {df['date'].iloc[-1].date()}, "
      f"mean {df['q'].mean():.1f} m3/s")
