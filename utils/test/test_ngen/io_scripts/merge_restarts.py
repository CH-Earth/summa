#!/usr/bin/env python3
"""
Merge NGEN restart files from multiple GRUs into a single SUMMA restart file with same GRU order as a given (old) restart file.
Usage:
  python3 merge_restarts.py init_cond.nc restart_pattern date_string
"""
import argparse
import xarray as xr
import os

def drop_hru_if_has_gru(ds, strategy="first"):
    """
    For every data variable that has both 'gru' and 'hru' dims,
    collapse the 'hru' dim according to strategy and replace the var.
    strategy: "first" | "mean" | "sum"
    """
    for vn in list(ds.data_vars):
        da = ds[vn]
        dims = tuple(str(d) for d in da.dims)
        if "gru" in dims and "hru" in dims:
            if strategy == "first":
                new = da.isel(hru=0)
            elif strategy == "mean":
                new = da.mean(dim="hru")
            elif strategy == "sum":
                new = da.sum(dim="hru")
            else:
                raise ValueError(f"unknown strategy {strategy}")
            ds[vn] = new
    # remove any leftover coordinate variables (optional)
    ds = ds.reset_coords(drop=True)
    return ds



def main():
    p = argparse.ArgumentParser()
    p.add_argument("init_cond", help="init_cond.nc (initial condition file with desired GRU order)")
    p.add_argument("restart_pattern", help="pattern for GRU restart files, e.g. gauge_01073000/simulations/run_1/SUMMA/run_1_restart_{}_G{}-{}.nc")
    p.add_argument("date_string", help="date string to identify restart files, e.g. 2015123023")
    args = p.parse_args()


    file_ds = xr.open_dataset(args.init_cond)
    hru_ids = file_ds['hruId'].values.tolist()

    # create empty list to hold datasets
    ds_list = []
    for i, hru_id in enumerate(hru_ids):
        restart_file = args.restart_pattern.format(args.date_string,str(i+1),str(i+1))
        if not os.path.isfile(restart_file):
            raise RuntimeError(f"Restart file {restart_file} for HRU ID {hru_id} not found.")
        ds = xr.open_dataset(restart_file)
        ds_list.append(ds)

    if not ds_list:
        raise RuntimeError("No restart files found.")

    # union of variable names across all datasets
    all_vars = sorted(set().union(*(set(ds.data_vars.keys()) for ds in ds_list)))

    out_vars = {}
    for var in all_vars:
        # choose concat dim for this variable
        choose_dim = None
        for ds in ds_list:
            if var in ds and "hru" in ds[var].dims:
                choose_dim = 'hru'
                break
        if choose_dim is None:
            for ds in ds_list:
                if var in ds and "gru" in ds[var].dims:
                    choose_dim = 'gru'
                    break

        # prepare per-file DataArray list, ensuring a length-1 concat dim and coord value from hru_ids
        arrs = []
        for idx, ds in enumerate(ds_list):
            this_id = hru_ids[idx]
            if var in ds:
                da = ds[var]
                dims = tuple(str(d) for d in da.dims)
                if choose_dim in dims:
                    da = da.expand_dims({}) if False else da  # no-op; keep shape
                    # ensure coord for the dim is the desired id
                    try:
                        da = da.assign_coords({choose_dim: (choose_dim, [this_id])})
                    except Exception:
                        # some variables may not allow assign_coords; construct new DataArray
                        da = xr.DataArray(da.values, dims=da.dims, coords={choose_dim: [this_id]}, attrs=da.attrs)
                else:
                    # add leading concat dim
                    da = da.expand_dims({choose_dim: [this_id]})
                    da = da.assign_coords({choose_dim: (choose_dim, [this_id])})
            else:
                raise RuntimeError(f"Variable {var} not found in restart file for HRU ID {this_id}")

            arrs.append(da)

        # concat for this variable
        try:
            conc = xr.concat(arrs, dim=choose_dim, combine_attrs='override')
        except Exception:
            # fallback: try concat by building numpy stack then DataArray
            conc = xr.concat(arrs, dim=choose_dim, combine_attrs='override')
        out_vars[var] = conc

    # assemble final dataset from per-variable DataArrays
    merged_ds = xr.Dataset(data_vars=out_vars)
    out_path = args.restart_pattern.replace("{}_G{}-{}.nc", f"merged_{args.date_string}.nc")
    merged_ds.to_netcdf(out_path)
    print(f"Wrote merged restart to {out_path}")

if __name__ == "__main__":
    main()