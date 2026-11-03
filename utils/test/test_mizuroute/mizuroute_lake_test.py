#!/usr/bin/env python3
"""Build, check and compare the SUMMA-mizuRoute lake test.

The bundled Provo domain has no open water, so `build` changes a copy of it:

    GRU A: land HRU -> reach 30, water HRU -> Doll lake, reach 31
    GRU B: land HRU -> reach 39, water HRU -> Doll lake, reach 40,
                                 water HRU -> HYPE lake, reach 41

All reaches form one chain, and a remap file gives each routing HRU the values
of its GRU. Each lake HRU is the lake only, and the land of its GRU drains to a
reach upstream of the lake.

Usage:
    mizuroute_lake_test.py build <settings_in> <work_dir>
    mizuroute_lake_test.py check <output.nc> <work_dir> <forcing.nc>
    mizuroute_lake_test.py standalone-input <output.nc> <work_dir>
    mizuroute_lake_test.py compare <output.nc> <standalone.nc>
"""

import os
import re
import shutil
import sys

import numpy as np
from netCDF4 import Dataset, date2num, num2date

# zero-based HRU indices in the bundled domain
GRU_A = {"land": 29, "water": [30]}
GRU_B = {"land": 38, "water": [39, 40]}
DOLL_LAKES = [30, 39]
HYPE_LAKES = [40]

# water land-cover class of each vegetation table
WATER_CLASS = {"USGS": 16, "MODIFIED_IGBP_MODIS_NOAH": 17}

SEGMENT_LENGTH_M = 5000.0
SEGMENT_SLOPE = 0.001

# multiconst.f90 latent heat of vaporization (J kg-1), density of water (kg m-3)
LH_VAP = 2501000.0
IDEN_WATER = 1000.0
# default Oudin parameters (deg C)
OUDIN_K1 = 100.0
OUDIN_K2 = 5.0

# tolerances: exact to rounding in double precision; the standalone history file is single precision
TOL_RATE = 1.0e-15  # m s-1
TOL_BALANCE = 1.0e-9  # fraction of the mean lake storage
TOL_STANDALONE = 1.0e-5  # fraction of the largest value at the reach


def copy_nc(src, dst, drop=None, edit=None):
    """Copy a NetCDF file, dropping indices along named dimensions and replacing named variables."""
    drop = drop or {}
    edit = edit or {}
    with Dataset(src) as a, Dataset(dst, "w", format="NETCDF4") as b:
        for name, dim in a.dimensions.items():
            b.createDimension(name, len(dim) - len(drop.get(name, [])))
        for name, var in a.variables.items():
            out = b.createVariable(name, var.datatype, var.dimensions)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            data = edit[name] if name in edit else var[:]
            for axis, dim in enumerate(var.dimensions):
                if dim in drop:
                    data = np.delete(data, drop[dim], axis=axis)
            out[:] = data


def build(settings_in, work_dir):
    settings = os.path.join(work_dir, "settings")
    mizu = os.path.join(work_dir, "mizu")
    os.makedirs(settings, exist_ok=True)
    os.makedirs(mizu, exist_ok=True)
    for name in os.listdir(settings_in):
        path = os.path.join(settings_in, name)
        if os.path.isfile(path):
            shutil.copy(path, settings)

    with open(os.path.join(settings_in, "modelDecisions.txt")) as f:
        table = re.search(r"^vegeParTbl\s+(\S+)", f.read(), re.M).group(1)
    water_class = WATER_CLASS[table]

    with Dataset(os.path.join(settings_in, "attributes.nc")) as a:
        gru_id = np.array(a["gruId"][:], dtype="i8")
        hru_id = np.array(a["hruId"][:], dtype="i8")
        hru2gru = np.array(a["hru2gruId"][:], dtype="i8")
        area = np.array(a["HRUarea"][:], dtype="f8")
        veg = np.array(a["vegTypeIndex"][:])

    # water HRUs join the GRU of their land HRU and those GRUs are dropped
    dropped = []
    for group in (GRU_A, GRU_B):
        for ix in group["water"]:
            veg[ix] = water_class
            dropped.append(int(np.where(gru_id == hru2gru[ix])[0][0]))
            hru2gru[ix] = hru2gru[group["land"]]

    copy_nc(os.path.join(settings_in, "attributes.nc"), os.path.join(settings, "attributes.nc"),
            drop={"gru": dropped}, edit={"vegTypeIndex": veg, "hru2gruId": hru2gru})
    copy_nc(os.path.join(settings_in, "trialParams.nc"), os.path.join(settings, "trialParams.nc"),
            drop={"gru": dropped})

    n = len(hru_id)
    seg_id = np.arange(1, n + 1, dtype="i8")
    is_lake = np.zeros(n, dtype="i4")
    lake_type = np.zeros(n, dtype="i4")
    is_lake[DOLL_LAKES + HYPE_LAKES] = 1
    lake_type[DOLL_LAKES] = 1
    lake_type[HYPE_LAKES] = 3

    with Dataset(os.path.join(mizu, "topology.nc"), "w", format="NETCDF4") as d:
        d.createDimension("hru", n)
        d.createDimension("seg", n)

        def put(name, dim, data, dtype):
            d.createVariable(name, dtype, (dim,))[:] = data

        put("hruId", "hru", hru_id, "i8")
        put("area", "hru", area, "f8")
        put("hruToSegId", "hru", seg_id, "i8")
        put("segId", "seg", seg_id, "i8")
        put("downSegId", "seg", np.concatenate([seg_id[1:], [0]]), "i8")
        put("length", "seg", np.full(n, SEGMENT_LENGTH_M), "f8")
        put("slope", "seg", np.full(n, SEGMENT_SLOPE), "f8")
        put("islake", "seg", is_lake, "i4")
        put("LakeModelType", "seg", lake_type, "i4")
        for name, value in (("D03_MaxStorage", 5.0e6), ("D03_Coefficient", 0.1), ("D03_Power", 1.5), ("D03_S0", 0.0),
                            ("HYP_E_emr", 10.0), ("HYP_E_lim", 8.0), ("HYP_E_min", 2.0), ("HYP_E_zero", 0.0),
                            ("HYP_Qrate_emr", 50.0), ("HYP_Erate_emr", 1.5), ("HYP_Qrate_prim", 100.0),
                            ("HYP_Qrate_amp", 0.3), ("HYP_Qrate_phs", 100.0), ("HYP_prim_F", 1.0),
                            ("HYP_Qsim_mode", 1.0)):
            put(name, "seg", np.full(n, value), "f8")
        put("HYP_A_avg", "seg", area, "f8")

    with Dataset(os.path.join(mizu, "remap.nc"), "w", format="NETCDF4") as d:
        d.createDimension("polyid", n)
        d.createDimension("data", n)
        d.createVariable("polyid", "i8", ("polyid",))[:] = hru_id
        d.createVariable("nOverlaps", "i4", ("polyid",))[:] = np.ones(n, dtype="i4")
        d.createVariable("qhru_id", "i8", ("data",))[:] = hru2gru
        d.createVariable("weight", "f8", ("data",))[:] = np.ones(n)

    print(f"built lake domain ({table}, water class {water_class}): {n} routing HRUs on {n - len(dropped)} GRUs,")
    print(f"  Doll lakes at reaches {seg_id[DOLL_LAKES]}, HYPE lake at reach {seg_id[HYPE_LAKES]}")


def read_output(output):
    with Dataset(output) as o:
        out = {
            "precip": np.array(o["lake_precip"][:], dtype="f8"),
            "evap": np.array(o["lake_evap"][:], dtype="f8"),
            "Q": np.array(o["Q_reach"][:, :, 0], dtype="f8"),
            "V": np.array(o["V_reach"][:, :, 0], dtype="f8"),
            "q_basin": np.array(o["q_basin"][:], dtype="f8"),
            "runoff": np.array(o["averageRoutedRunoff"][:], dtype="f8"),
            "seg": np.array(o["seg"][:]),
            "dates": num2date(o["time"][:], o["time"].units),
        }
    return out


def oudin(sw, airtemp):
    """Oudin potential evaporation (m s-1), written out separately from the Fortran."""
    pe = (sw / (LH_VAP * IDEN_WATER)) * ((airtemp - 273.15 + OUDIN_K2) / OUDIN_K1)
    return np.maximum(pe, 0.0)


def check(output, work_dir, forcing):
    out = read_output(output)
    fail = False

    def report(ok, text):
        nonlocal fail
        fail = fail or not ok
        print(("  ok    " if ok else "  FAIL  ") + text)

    with Dataset(os.path.join(work_dir, "settings", "attributes.nc")) as a:
        gru_id = np.array(a["gruId"][:], dtype="i8")
        hru_id = np.array(a["hruId"][:], dtype="i8")
        hru2gru = np.array(a["hru2gruId"][:], dtype="i8")
        area = np.array(a["HRUarea"][:], dtype="f8")
    with Dataset(forcing) as f:
        hours = np.rint(date2num(out["dates"], f["time"].units)).astype(int)
        ix_time = [int(np.where(np.array(f["time"][:]).astype(int) == h)[0][0]) for h in hours]
        hru_file = np.array(f["hruId"][:], dtype="i8")
        ix_hru = {int(h): int(np.where(hru_file == h)[0][0]) for h in hru_id}

        def forc(name, ix):
            return np.array(f[name][ix_time, ix_hru[int(hru_id[ix])]], dtype="f8")

        gru_of = lambda ix: int(np.where(gru_id == hru2gru[ix])[0][0])
        report(np.array_equal(out["q_basin"], out["runoff"][:, [gru_of(i) for i in range(len(hru_id))]]),
               "every routing HRU receives its GRU's runoff through the remap")

        # lake precipitation and evaporation, area-weighted across the water HRUs of the GRU
        water_grus = set()
        for group in (GRU_A, GRU_B):
            ixs = group["water"]
            g = gru_of(group["land"])
            water_grus.add(g)
            total = area[ixs].sum()
            p = sum(forc("pptrate", i) * area[i] for i in ixs) / total / IDEN_WATER
            e = sum(oudin(forc("SWRadAtm", i), forc("airtemp", i)) * area[i] for i in ixs) / total
            dp, de = np.abs(out["precip"][:, g] - p).max(), np.abs(out["evap"][:, g] - e).max()
            report(dp < TOL_RATE and de < TOL_RATE and p.max() > 0.0 and e.max() > 0.0,
                   f"GRU {gru_id[g]} ({len(ixs)} water HRU): max|P error| {dp:.2e}, max|E error| {de:.2e} m/s "
                   f"(peaks {p.max():.2e}, {e.max():.2e})")
    others = [g for g in range(len(gru_id)) if g not in water_grus]
    report(np.abs(out["precip"][:, others]).max() == 0.0 and np.abs(out["evap"][:, others]).max() == 0.0,
           "GRUs without a water HRU supply no lake precipitation or evaporation")

    # water balance of each lake, and flow is smoothed compared with its inflow
    dt = 3600.0
    seg = list(out["seg"])
    for ix in DOLL_LAKES + HYPE_LAKES:
        i = seg.index(ix + 1)
        up = seg.index(ix)
        g = gru_of(ix)
        to_volume = area[ix] * dt
        inflow = (out["Q"][1:, up] * dt
                  + (out["q_basin"][1:, ix] + out["precip"][1:, g] - out["evap"][1:, g]) * to_volume
                  - out["Q"][1:, i] * dt)
        residual = np.abs(np.diff(out["V"][:, i]) - inflow).max()
        rel = residual / out["V"][:, i].mean()
        kind = "HYPE" if ix in HYPE_LAKES else "Doll"
        report(rel < TOL_BALANCE, f"{kind} lake at reach {ix + 1}: water balance closes to {residual:.2e} m3 ({rel:.1e} of storage)")
        report(out["Q"][:, i].std() < out["Q"][:, up].std(),
               f"{kind} lake at reach {ix + 1}: outflow is smoother than inflow "
               f"(std {out['Q'][:, i].std():.1f} vs {out['Q'][:, up].std():.1f} m3/s)")

    if np.isnan(out["Q"]).any() or np.isnan(out["V"]).any() or (out["Q"] < 0.0).any():
        report(False, "NaN or negative streamflow")
    return 1 if fail else 0


def standalone_input(output, work_dir):
    """Write what the coupled run sent to mizuRoute as a standalone mizuRoute forcing file."""
    out = read_output(output)
    with Dataset(os.path.join(work_dir, "settings", "attributes.nc")) as a:
        gru_id = np.array(a["gruId"][:], dtype="i8")
    units = "hours since 2017-01-01 00:00:00"
    path = os.path.join(work_dir, "standalone", "forcing.nc")
    with Dataset(path, "w", format="NETCDF4") as d:
        d.createDimension("time", None)
        d.createDimension("hru", len(gru_id))
        t = d.createVariable("time", "f8", ("time",))
        t.units = units
        t.calendar = "standard"
        t[:] = np.rint(date2num(out["dates"], units))
        d.createVariable("hruId", "i8", ("hru",))[:] = gru_id
        for name, data in (("runoff", out["runoff"]), ("precip", out["precip"]), ("evap", out["evap"])):
            var = d.createVariable(name, "f8", ("time", "hru"), fill_value=-9999.0)
            var.units = "m/s"
            var[:] = data
    print(f"wrote {path}")


def compare(output, standalone):
    out = read_output(output)
    with Dataset(standalone) as s:
        reach = list(np.array(s["reachID"][:]))
        q = np.array(s["KWroutedRunoff"][:], dtype="f8")
        v = np.array(s["KWvolume"][:], dtype="f8")
    order = [reach.index(r) for r in out["seg"]]
    q, v = q[:, order], v[:, order]
    if q.shape != out["Q"].shape:
        print(f"  FAIL  time steps or reaches differ: {q.shape} standalone, {out['Q'].shape} coupled")
        return 1

    fail = False
    for name, coupled, alone in (("streamflow", out["Q"], q), ("storage", out["V"], v)):
        err = np.abs(coupled - alone).max(axis=0)
        rel = err / np.maximum(np.abs(coupled).max(axis=0), np.finfo("f8").tiny)
        worst = int(np.argmax(rel))
        ok = rel.max() < TOL_STANDALONE
        fail = fail or not ok
        print(("  ok    " if ok else "  FAIL  ")
              + f"{name}: largest difference {rel.max():.2e} of the reach maximum, at reach {out['seg'][worst]}")
    return 1 if fail else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    command = sys.argv[1]
    if command == "build":
        build(sys.argv[2], sys.argv[3])
    elif command == "check":
        sys.exit(check(sys.argv[2], sys.argv[3], sys.argv[4]))
    elif command == "standalone-input":
        standalone_input(sys.argv[2], sys.argv[3])
    elif command == "compare":
        sys.exit(compare(sys.argv[2], sys.argv[3]))
    else:
        sys.exit(__doc__)
