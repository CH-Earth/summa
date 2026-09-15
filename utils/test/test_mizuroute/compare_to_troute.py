#!/usr/bin/env python3
"""Compare SUMMA-mizuRoute routed flow to t-route routed flow on the same network.

Prerequisite: run both routing schemes over the same real hydrofabric domain
and the same simulated period:

  - mizuRoute: test_mizuroute_provo_real_network.sh, which converts the
    hydrofabric with hydrofabric_to_topology.py and writes Q_reach.
  - t-route: ../test_ngen/provo_run.sh (run from the ngen repo root -- this
    needs ngen built and a python env with nwm_routing importable, see
    ../test_ngen/readme.md), which writes troute_output_*.nc into
    domain_provo/simulations/.

Both routing schemes key their output by the hydrofabric's numeric catchment
id (mizuRoute's segId, set by hydrofabric_to_topology.py, equals t-route's
feature_id), so no id crosswalk is needed -- this script matches directly on
that id.

This is not expected to match closely: mizuRoute is run here with the
kinematic wave method on reach length/slope alone, while t-route defaults to
a Muskingum-Cunge-like scheme using real channel geometry (width, Manning's
n) that the SUMMA-mizuRoute coupling does not currently pass through. Use
this to sanity-check that both are routing the same real network (same
reaches respond, similar timing and volumes), not to validate mizuRoute
against t-route as a reference.

Usage:
    compare_to_troute.py <mizuroute_timestep.nc> <troute_output_dir_or_file> [topology.nc]
"""

import glob
import os
import sys

import numpy as np
from netCDF4 import Dataset, num2date


def load_mizuroute(path):
    with Dataset(path) as d:
        seg = np.array(d.variables["seg"][:], dtype="i8")
        q = np.ma.filled(d.variables["Q_reach"][:, :, 0].astype("f8"), np.nan)  # (time, seg)
        time_var = d.variables["time"]
        time = num2date(time_var[:], time_var.units, only_use_cftime_datetimes=False)
    # num2date introduces sub-second float jitter; round to the nearest second
    return seg, np.array([np.datetime64(t, "s") for t in time]), q


def load_troute(path):
    files = sorted(glob.glob(os.path.join(path, "troute_output_*.nc"))) if os.path.isdir(path) else [path]
    if not files:
        sys.exit(f"no troute_output_*.nc files found in {path}")

    feature_id = None
    times, flows = [], []
    for f in files:
        with Dataset(f) as d:
            fid = np.array(d.variables["feature_id"][:], dtype="i8")
            if feature_id is None:
                feature_id = fid
            elif not np.array_equal(feature_id, fid):
                sys.exit(f"feature_id order differs between t-route output files (offender: {f})")

            flow = np.ma.filled(d.variables["flow"][:].astype("f8"), np.nan)  # (feature_id, time)
            time_var = d.variables["time"]
            t0 = np.datetime64(time_var.units.replace("seconds since ", ""))
            t = t0 + np.array(time_var[:], dtype="int64").astype("timedelta64[s]")

        times.append(t)
        flows.append(flow)

    time = np.concatenate(times)
    flow = np.concatenate(flows, axis=1)  # (feature_id, time)
    order = np.argsort(time)
    return feature_id, time[order], flow[:, order]


def nse(sim, obs):
    mask = np.isfinite(sim) & np.isfinite(obs)
    if mask.sum() < 2:
        return np.nan
    sim, obs = sim[mask], obs[mask]
    denom = np.sum((obs - obs.mean()) ** 2)
    if denom == 0:
        return np.nan
    return 1.0 - np.sum((sim - obs) ** 2) / denom


def main(mizu_file, troute_path, topology_file=None):
    mizu_seg, mizu_time, mizu_q = load_mizuroute(mizu_file)
    troute_id, troute_time, troute_q = load_troute(troute_path)

    print(f"mizuRoute: {len(mizu_seg)} reaches, {len(mizu_time)} steps, "
          f"{mizu_time[0]} to {mizu_time[-1]}")
    print(f"t-route:   {len(troute_id)} reaches, {len(troute_time)} steps, "
          f"{troute_time[0]} to {troute_time[-1]}")

    common_ids = np.intersect1d(mizu_seg, troute_id)
    if len(common_ids) == 0:
        sys.exit("no reach ids in common between mizuRoute and t-route output")
    print(f"reaches in common: {len(common_ids)} / {len(mizu_seg)}")

    t_start = max(mizu_time[0], troute_time[0])
    t_end = min(mizu_time[-1], troute_time[-1])
    if t_start >= t_end:
        sys.exit("mizuRoute and t-route outputs do not overlap in time")

    # align both series onto mizuRoute's timestamps (both are written hourly here)
    common_time = mizu_time[(mizu_time >= t_start) & (mizu_time <= t_end)]
    mizu_idx = np.searchsorted(mizu_time, common_time)
    troute_idx = np.searchsorted(troute_time, common_time)
    troute_idx = np.clip(troute_idx, 0, len(troute_time) - 1)

    outlet_id = None
    if topology_file:
        with Dataset(topology_file) as t:
            seg_id = np.array(t.variables["segId"][:], dtype="i8")
            down_seg = np.array(t.variables["downSegId"][:], dtype="i8")
        outlet_id = int(seg_id[down_seg == 0][0])

    print()
    print(f"{'reach id':>12} {'mizu mean':>12} {'troute mean':>12} "
          f"{'mizu peak':>12} {'troute peak':>12} {'NSE':>8}")
    print("-" * 76)

    nse_values = []
    for rid in sorted(common_ids):
        m_series = mizu_q[mizu_idx, np.where(mizu_seg == rid)[0][0]]
        t_series = troute_q[np.where(troute_id == rid)[0][0], troute_idx]

        score = nse(m_series, t_series)
        nse_values.append(score)

        marker = " (outlet)" if rid == outlet_id else ""
        print(f"{rid:>12} {np.nanmean(m_series):>12.3f} {np.nanmean(t_series):>12.3f} "
              f"{np.nanmax(m_series):>12.3f} {np.nanmax(t_series):>12.3f} "
              f"{score:>8.3f}{marker}")

    nse_values = np.array(nse_values, dtype="f8")
    print()
    print(f"median NSE across {len(common_ids)} reaches: {np.nanmedian(nse_values):.3f}")
    print("(NSE of 1 is a perfect match; this compares two different routing")
    print(" schemes, so a modest score is expected -- see the module docstring)")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    main(*sys.argv[1:])
