#!/usr/bin/env python3
"""
Convert SUMMA forcing files to NGen forcing format or vice versa.
Usage:
  python3 convert_forcing.py input_file.nc output_file.nc typeOfConversion attributes_file.nc(optional)
"""
import argparse
import xarray as xr
import numpy as np
import re

def _parse_id(fid):
    # try direct int, else extract trailing digits, else first digits, else return original string
    try:
        return int(fid)
    except Exception:
        s = str(fid).strip()
        m = re.search(r'(\d+)$', s) or re.search(r'(\d+)', s)
        return int(m.group(1)) if m else s

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_file", help="input_file.nc (SUMMA or NGen forcing file)")
    p.add_argument("output_file", help="output_file.nc (converted forcing file)")
    p.add_argument("typeOfConversion", help="type of conversion: summa_to_ngen or ngen_to_summa")
    p.add_argument("attributes_file", nargs='?', help="attributes_file.nc (optional, for lat/lon info)", default=None)
    args = p.parse_args()

    # open input file
    f_ds = xr.open_dataset(args.input_file)

    # do conversion
    if args.typeOfConversion == "ngen_to_summa":
        f_var = "ids"
        if f_var not in f_ds:
            raise RuntimeError(f"NGen forcing file missing '{f_var}' variable")
        f_ids_raw = np.array([str(x) for x in f_ds[f_var].values])
        f_ids = np.array([_parse_id(fid) for fid in f_ids_raw])

        # make hru dimension name consistent with SUMMA and replace ids with parsed integer ids
        f_ds = f_ds.rename_dims({"catchment-id": "hru"})
        f_ds[f_var] = (("hru",), f_ids)

        # rename other variables to match SUMMA conventions
        f_ds = f_ds.rename_vars({
            "ids": "hruId",
            "precip_rate": "pptrate",
            "TMP_2maboveground": "airtemp",
            "SPFH_2maboveground": "spechum",
            "PRES_surface": "airpres",
            "DSWRF_surface": "SWRadAtm",
            "DLWRF_surface": "LWRadAtm",
        })

        # make time a coordinate instead of variable and set values to Time values as datetime64[ns]
        time_da = f_ds['Time']
        if time_da.ndim == 1:
            tvals = time_da.values.astype('datetime64[ns]')
        elif time_da.ndim == 2:
            # if Time is repeated per hru, take the first hru column (assumes identical)
            other_dims = [d for d in time_da.dims if d != 'time']
            if other_dims:
                tvals = time_da.isel({other_dims[0]: 0}).values.astype('datetime64[ns]')
            else:
                # fallback: assume shape (time, N) and take first column
                tvals = time_da.values[:, 0].astype('datetime64[ns]')
        else:
            raise RuntimeError(f"Unsupported Time ndim: {time_da.ndim}")
        f_ds = f_ds.assign_coords(time=('time', tvals))
        f_ds = f_ds.drop_vars('Time')

        # make windspd UV components into speed
        u = f_ds["UGRD_10maboveground"]
        v = f_ds["VGRD_10maboveground"]
        f_ds["windspd"] = np.sqrt(u**2 + v**2)
        f_ds = f_ds.drop_vars("UGRD_10maboveground")
        f_ds = f_ds.drop_vars("VGRD_10maboveground")

        # add data_step variable
        f_ds["data_step"] = ((), np.int64(3600))

        # add lat/lon variables as read from attributes (if attributes_file provided)
        if args.attributes_file:
            f_attr = xr.open_dataset(args.attributes_file)
            f_ds["latitude"] = f_attr["latitude"]
            f_ds["longitude"] = f_attr["longitude"]

        # ensure data variables use dims ('time','hru') not ('hru','time')
        for name, var in list(f_ds.data_vars.items()):
            dims = tuple(var.dims)
            if dims == ("hru", "time"):
                f_ds[name] = var.transpose("time", "hru")
            elif set(dims) >= {"time", "hru"} and dims.index("hru") < dims.index("time"):
                # generic reorder if hru appears before time
                new_order = tuple(d for d in dims if d != "hru") + ("hru",) if "time" not in dims else ("time","hru")
                try:
                    f_ds[name] = var.transpose(*new_order)
                except Exception:
                    pass

        # write output file
        f_ds.to_netcdf(args.output_file)
        print(f"Wrote converted SUMMA forcing file to {args.output_file}")

    elif args.typeOfConversion == "summa_to_ngen":
        f_var = "hruId"
        if f_var not in f_ds:
            raise RuntimeError(f"SUMMA forcing file missing '{f_var}' variable")
        f_ids_raw = f_ds[f_var].values
        f_ids = np.array([f"cat-{str(fid)}" for fid in f_ids_raw])

        # make hru dimension name consistent with NGen and replace hruId with string ids
        f_ds = f_ds.rename_dims({"hru": "catchment-id"})
        f_ds[f_var] = (("catchment-id"), f_ids)

        # rename other variables to match NGen conventions
        f_ds = f_ds.rename_vars({
            "hruId": "ids",
            "pptrate": "precip_rate",
            "airtemp": "TMP_2maboveground",
            "spechum": "SPFH_2maboveground",
            "airpres": "PRES_surface",
            "SWRadAtm": "DSWRF_surface",
            "LWRadAtm": "DLWRF_surface",
        })

        # make Time a variable from time coordinate tiled by catchment and drop time coordinate
        time_vals = f_ds.coords["time"].values.astype("datetime64[ns]")
        time_vals = np.asarray(time_vals, dtype="datetime64[ns]")
        ntime = time_vals.shape[0]
        n = f_ds.sizes.get("catchment-id", f_ds.sizes.get("hru", 1))

        # build 2-D Time in nanoseconds as float64 to match other files
        time_ns = time_vals.astype("datetime64[ns]").astype("int64").astype(np.float64)
        tiled = np.broadcast_to(time_ns[None, :], (n, ntime))  # shape (catchment-id, time)
        f_ds["Time"] = (("catchment-id", "time"), tiled)
        f_ds["Time"].attrs["units"] = "ns"
        # clear problematic encoding/attrs that can trigger packing/valid_range errors
        f_ds["Time"].encoding.clear()
        for bad in ("valid_range", "_FillValue", "missing_value", "scale_factor", "add_offset"):
            f_ds["Time"].attrs.pop(bad, None)

        # make windspd into UV components, all north wind for now
        windspd = f_ds["windspd"]
        f_ds["UGRD_10maboveground"] = windspd
        f_ds["VGRD_10maboveground"] = xr.zeros_like(windspd)
        f_ds = f_ds.drop_vars("windspd")

        # remove SUMMA-only variables (if present)
        for maybe in ("latitude", "longitude", "data_step"):
                   if maybe in f_ds:
                        f_ds = f_ds.drop_vars(maybe)

        # ensure data variables use dims ('catchment-id','time') not ('time','catchment-id')
        for name, var in list(f_ds.data_vars.items()):
            dims = tuple(var.dims)
            if dims == ("time", "catchment-id"):
                f_ds[name] = var.transpose("catchment-id", "time")
            elif set(dims) >= {"time", "catchment-id"} and dims.index("time") < dims.index("catchment-id"):
                # generic reorder if time appears before catchment-id
                new_order = tuple(d for d in dims if d != "time") + ("time",) if "catchment-id" not in dims else ("catchment-id","time")
                try:
                    f_ds[name] = var.transpose(*new_order)
                except Exception:
                    pass

        # write output file
        f_ds.to_netcdf(args.output_file)
        print(f"Wrote NGen forcing to {args.output_file}")
            
    else:
        raise RuntimeError("typeOfConversion must be either 'summa_to_ngen' or 'ngen_to_summa'")






if __name__ == "__main__":
    main()