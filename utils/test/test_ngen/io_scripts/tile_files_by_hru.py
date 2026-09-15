#!/usr/bin/env python3
"""
Tile/reindex SUMMA file with one HRU by hruId to match ids/order in NGEN forcing.nc.
Usage:
  python3 tile_files_by_hru.py old_file.nc forcing.nc typeOfInput folder_to_other_attribs (optional) 
"""
import argparse
import xarray as xr
import numpy as np
import re
import os
from ast import literal_eval

def _parse_id(fid):
    # try direct int, else extract trailing digits, else first digits, else return original string
    try:
        return int(fid)
    except Exception:
        s = str(fid).strip()
        m = re.search(r'(\d+)$', s) or re.search(r'(\d+)', s)
        return int(m.group(1)) if m else s
    
def _parse_namelist(text):
    groups = {}
    for m in re.finditer(r'&(\w+)(.*?)/', text, re.S):
        name, body = m.group(1), m.group(2)
        vals = {}
        for line in re.split(r'[\r\n]+', body):
            line = re.sub(r'!.*$', '', line).strip()   # drop comments
            if not line: 
                continue
            for kv in line.split(','):
                if '=' not in kv: 
                    continue
                k, v = map(str.strip, kv.split('=', 1))
                # try to convert to number/list/string
                try:
                    parsed = literal_eval(v)
                except Exception:
                    parsed = v.strip().strip('"').strip("'")
                vals[k] = parsed
        groups[name] = vals
    return groups

def _parse_kv_text(text):
    from ast import literal_eval
    import re

    out = {}
    for line in text.splitlines():
        line = re.sub(r'#.*$', '', line).strip()            # drop comments starting with #
        line = re.sub(r'\[.*$', '', line).strip()           # drop units starting with [
        if not line: 
            continue
        if '=' not in line:
            continue
        k, v = map(str.strip, line.split('=', 1))
        # normalize booleans used in config files
        vl = v.strip().strip('"').strip("'")
        if vl.upper() in ("TRUE", "T", "YES", "Y", "ON"):
            out[k] = True
            continue
        if vl.upper() in ("FALSE", "F", "NO", "N", "OFF"):
            out[k] = False
            continue
        # try numeric / list parsing
        try:
            out[k] = literal_eval(vl)
        except Exception:
            # fallback: try int/float then raw string
            try:
                out[k] = int(vl)
            except Exception:
                try:
                    out[k] = float(vl)
                except Exception:
                    out[k] = vl
    return out

def _ensure_hru_var(ds, name, n, fill=np.nan, dtype=float):
    if name not in ds:
        arr = np.full((n,), fill, dtype=dtype)
        ds[name] = xr.DataArray(arr, dims=("hru",), attrs={})
    return ds

def _fixSoilVar(ds,var_name):
    """Fix soil variable to have midSoil dimension instead of midToto."""
    da = ds[var_name]

    if 'midToto' in da.dims:
        # rename the existing dimension
        ds['mLayerMatricHead'] = da.rename({'midToto': 'midSoil'})
    else:
        # create a new midSoil dim of the same length as midToto (fallback to 1 if missing)
        size = ds.dims.get('midToto', 1)
        ds[var_name] = da.expand_dims({'midSoil': size})
    return ds

def _fixWaterContentVar(ds):
    """Fix water content variable to be small enough."""
    # if mLayerVolFracLiq is close to default value then lower it
    default_water_content = 0.45 # assumed default
    ds['mLayerVolFracLiq'].values[ds['mLayerVolFracLiq'].values > default_water_content-0.] = default_water_content - 0.15
    return ds

def main():
    p = argparse.ArgumentParser()
    p.add_argument("old_file", help="old_file.nc (with single HRU)")
    p.add_argument("forcing", help="forcing.nc (contains target HRU ids/order)")
    p.add_argument("typeOfInput", help="type of input files in other_folder: param or attrib ")
    p.add_argument("other_folder", nargs='?', default=None, help="file folder containing other model old_file to search (optional)")
    args = p.parse_args()

    file_ds = xr.open_dataset(args.old_file)
    f_ds = xr.open_dataset(args.forcing)
    typeOfInput = args.typeOfInput.lower()
    if typeOfInput not in ("param", "attrib", "init", "force"):
        raise RuntimeError(f"typeOfInput must be 'param' or 'attrib' or 'init' or 'force', got '{typeOfInput}'")
    
    if typeOfInput == "init":
        file_ds = _fixSoilVar(file_ds, 'mLayerMatricHead')
        file_ds = _fixWaterContentVar(file_ds)
    
    out_path = args.old_file.replace(".nc", f"_tiled_by_hru.nc")

    # forcing ids (use as labels for the new hru coordinate)
    f_var = "ids"
    if f_var not in f_ds:
        raise RuntimeError(f"forcing file missing '{f_var}' variable")
    f_ids_raw = np.array([str(x) for x in f_ds[f_var].values])
    f_ids = np.array([_parse_id(fid) for fid in f_ids_raw])
    n = f_ids.size

    # tile old_file to match forcing ids/order
    data_vars = {}
    for vn, var in file_ds.data_vars.items():
        var_dims = tuple(str(d) for d in var.dims)
        vals = var.values

        # prefer expanding an existing 'hru' axis, otherwise 'gru', otherwise add leading 'hru'
        if "hru" in var_dims:
            axis = var_dims.index("hru")
            size = var.sizes.get("hru", 1)
            if size == n:
                arr_vals = vals
            else:
                arr_vals = np.repeat(vals, n, axis=axis)
            arr = xr.DataArray(arr_vals, dims=var_dims, attrs=var.attrs)

        elif "gru" in var_dims:
            axis = var_dims.index("gru")
            size = var.sizes.get("gru", 1)
            if size == n:
                arr_vals = vals
            else:
                arr_vals = np.repeat(vals, n, axis=axis)
            arr = xr.DataArray(arr_vals, dims=var_dims, attrs=var.attrs)

        else:
            # keep same, no additional dimensions
            arr = xr.DataArray(var.values, dims=var_dims, attrs=var.attrs)

        data_vars[vn] = arr

    # preserve coordinate variables from original file_ds
    coords = {}
    for k, v in file_ds.coords.items():
        coords[k] = (tuple(str(d) for d in v.dims), v.values)

    tiled = xr.Dataset(data_vars=data_vars, coords=coords, attrs=file_ds.attrs)

    # remove all coordinate variables except 'time' (if any)
    to_drop = [c for c in list(tiled.coords) if c != 'time']
    if to_drop:
        tiled = tiled.drop_vars(to_drop)

    # make hruId and gruId be ids from forcing file
    if "hruId" in tiled:
        tiled['hruId'].values = f_ids
    if "gruId" in tiled:
        tiled['gruId'].values = f_ids

    # redo hru2gru mapping if present
    if "hru2gruId" in tiled:
        tiled['hru2gruId'].values = f_ids

    # set old_file from distributed files, all text files with matching hruId names
    NOAH_folder = os.path.join(args.other_folder, "NOAH")
    PET_folder = os.path.join(args.other_folder, "PET")
    CFE_folder = os.path.join(args.other_folder, "CFE")
            
    # default NOAH decisions
    #  precip_phase_option               = 1 # Jordan (1991) SNTHERM equation
    #  snow_albedo_option                = 1 # 1 is BATS==varDecay
    #  dynamic_veg_option                = 4 # use table LAI; use maximum vegetation fraction
    #  runoff_option                     = 3 # original surface and subsurface runoff (free drainage)
    #  drainage_option                   = 8 # original subsurface runoff (free drainage) with dynamic VIC runoff
    #  frozen_soil_option                = 1 # linear effects, more permeable
    #  dynamic_vic_option                = 1 # 1 is Philip scheme, 2 is Green-Ampt
    #  radiative_transfer_option         = 3 # modified two-stream (gap = F(solar angle, 3D structure ...)<1-FVEG)
    #  sfc_drag_coeff_option             = 1 # surface layer drag coeff M-O theory
    #  canopy_stom_resist_option         = 1 # Ball-Berry
    #  snowsoil_temp_time_option         = 3 # semi-implicit; flux top boundary condition, FSNO for TS calculation
    #  soil_temp_boundary_option         = 2 # zero flux bottom boundary condition
    #  supercooled_water_option          = 1 # no iteration for supercooled water
    #  stomatal_resistance_option        = 1 # NOAH, thresholded linear function of volumetric liquid water content
    #  evap_srfc_resistance_option       = 4 # Sakaguchi and Zeng (2009) for bare soil evaporation resistance, rsurf = rsurf_snow for snow
    #  subsurface_option                 = 2 # one-way coupled hydrostatic
    for fn in os.listdir(NOAH_folder):
        file_path = os.path.join(NOAH_folder, fn)
        base_name, ext = os.path.splitext(fn)
        # id should be number matching hruIds
        try:
            fid = _parse_id(base_name)
        except Exception:
            continue
        if fid in f_ids:
            idx = np.where(f_ids == fid)[0][0]
            nml = _parse_namelist(open(file_path, 'r').read())
            if typeOfInput == "param":
                # set parameters
                tiled = _ensure_hru_var(tiled, 'tempCritRain', n, fill=np.nan)
                tiled['tempCritRain'].values[idx] = nml['forcing']['rain_snow_thresh']+273.15  # convert C to K
            elif typeOfInput == "attrib":
                # NOTE: ideally would have HRUarea, but since GRU and HRU are 1-1, HRUarea is not needed in summa code
                tiled['latitude'].values[idx] = nml['location']['lat']
                tiled['longitude'].values[idx] = nml['location']['lon']
                #tiled['aspect'].values[idx] = nml['location']['azimuth'] # seems to be always 0
                tiled['mHeight'].values[idx] = nml['forcing']['ZREF'] # reference height for wind
                tiled['vegTypeIndex'].values[idx] = nml['structure']['vegtyp']
                tiled['soilTypeIndex'].values[idx] = nml['structure']['isltyp']
                #tiled['tan_slope'].values[idx] = np.tan(np.deg2rad(nml['location']['terrain_slope'])) # compute tan_slope from slope in degrees, seems to be always 0
            elif typeOfInput == "force":
                tiled['latitude'].values[idx] = nml['location']['lat']
                tiled['longitude'].values[idx] = nml['location']['lon']

    # PET 
    # uses logBelowCanopy wind decision
    for fn in os.listdir(PET_folder):
        file_path = os.path.join(PET_folder, fn)
        base_name, ext = os.path.splitext(fn)
        # id should be number matching hruIds
        try:
            fid = _parse_id(base_name)
        except Exception:
            continue
        if fid in f_ids:
            idx = np.where(f_ids == fid)[0][0]
            kvs = _parse_kv_text(open(file_path, 'r').read())
            if typeOfInput == "param":
                tiled = _ensure_hru_var(tiled, 'heightCanopyTop', n, fill=np.nan)
                tiled['heightCanopyTop'].values[idx] = kvs['vegetation_height_m']
                #tiled['zpdFraction'].values[idx] = kvs['zero_plane_displacement_height_m']/kvs['vegetation_height'] # seems too small
            elif typeOfInput == "attrib":
                tiled['elevation'].values[idx] = kvs['site_elevation_m']

    # CFE
    for fn in os.listdir(CFE_folder):
        file_path = os.path.join(CFE_folder, fn)
        base_name, ext = os.path.splitext(fn)
        # id should be number matching hruIds
        try:
            fid = _parse_id(base_name)
        except Exception:
            continue
        if fid in f_ids:
            idx = np.where(f_ids == fid)[0][0]
            kvs = _parse_kv_text(open(file_path, 'r').read())
            if typeOfInput == "param":
                continue # no parameters to set from CFE files, these available but contradict NOAH table values
                #tiled = _ensure_hru_var(tiled, 'theta_sat', n, fill=np.nan)
                #tiled = _ensure_hru_var(tiled, 'k_soil', n, fill=np.nan)
                #tiled = _ensure_hru_var(tiled, 'critSoilWilting', n, fill=np.nan)
                #tiled['theta_sat'].values[idx] = kvs['soil_params.smcmax'] # effective porosity [V/V]
                #tiled['k_soil'].values[idx] = kvs['soil_params.satdk'] # saturated hydraulic conductivity [m s-1]
                #tiled['critSoilWilting'].values[idx] = kvs['soil_params.wltsmc'] #  wilting point soil moisture content [V/V]

    if typeOfInput == "init":
        tiled['dt_init'].values = np.full(tiled['dt_init'].shape, 3600.0) # default 1 hour time step for initialization files

    # write result
    tiled.to_netcdf(out_path)
    print(f"Wrote tiled file to {out_path}")


if __name__ == "__main__":
    main()