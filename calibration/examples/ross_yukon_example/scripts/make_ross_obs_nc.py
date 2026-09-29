"""Convert the Ross daily streamflow CSV into the netCDF format J() reads
(same layout as CAN_05BB001_daily_flow_observations.nc)."""
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

obs_dir = Path(__file__).resolve().parent.parent / "observations"
df = pd.read_csv(obs_dir / "Ross_new_test2_streamflow_processed.csv")
df["datetime"] = pd.to_datetime(df["datetime"], format="%m/%d/%Y %H:%M")
q = df.set_index("datetime")["discharge_cms"].astype(float).sort_index()
q = q[~q.index.duplicated()].resample("D").mean()   # one value per day, gaps -> NaN

minutes = ((q.index - pd.Timestamp("1950-01-01")) // pd.Timedelta(minutes=1)).astype("int64")
ds = xr.Dataset(
    {"q_obs": ("time", q.values, {"units": "m3 s-1", "long_name": "observed streamflow values"})},
    coords={"time": ("time", minutes.values, {"standard_name": "time",
            "units": "minutes since 1950-01-01", "calendar": "proleptic_gregorian"})})
out = obs_dir / "Ross_daily_flow_observations.nc"
ds.to_netcdf(out, format="NETCDF4",
             encoding={"q_obs": {"_FillValue": np.nan}, "time": {"dtype": "int64"}})
print(f"wrote {out}\n  {q.index[0].date()} to {q.index[-1].date()}, {len(q)} days, "
      f"{int(q.isna().sum())} missing, {int((q < 0).sum())} negative")
