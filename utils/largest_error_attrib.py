#!/usr/bin/env python

# python largest_error_attrib.py [method_name] [stat]

import numpy as np
import xarray as xr
from pathlib import Path

top_fold    = '/home/avanb/scratch/'
attr_fold = '/home/avanb/TestScripts/settings/'
nBig = 10
do_rel = True # plot relative to the benchmark simulation

testing = False
if testing:
    top_fold = '/Users/amedin/Research/USask/test_py/'
    attr_fold = '/Users/amedin/Research/USask/test_py/settings/'
    method_name= 'be1'
    stat = 'rmse'
else:
    import sys
    method_name = sys.argv[1]
    stat = sys.argv[2]

des_dir =  top_fold + 'statistics'
des_dir = Path(des_dir)

settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name + '_hrly_diff_stats_{}.nc'
viz_fil = viz_fil.format(','.join(settings))
src_file =  des_dir / viz_fil
plot_vars = settings.copy()
attr_fil = Path(attr_fold) / 'attributes.nc'

# Open the netCDF file with RMSE data
summa = xr.open_dataset(src_file)

for var in plot_vars:
    stat0 = stat
    if stat == 'rmse' or stat == 'kgem': 
        if var == 'wallClockTime': stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz':
        if var == 'wallClockTime': stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe': 
        if var == 'wallClockTime': stat0 = 'amax'
        statr = 'amax_ben'

    # Get the variable from the netCDF file
    s = summa[var].sel(stat=stat0)
    if do_rel: 
        s_rel = summa[var].sel(stat=statr)
        if var != 'wallClockTime': s = s/s_rel

    # Mask the finite values of the variable
    mask = np.isfinite(s)
    s = s[mask]

    # Get the indices of the largest nBig values
    big_indices = s.argsort()[-nBig:]
    
    # Get the largest nBig values
    val_big = s[big_indices.values]
    
    # Get the hru coordinate of the largest nBig values
    hru_big = s[big_indices.values].hru.values
    
    # Print all the values of the biggest rmse hru
    if do_rel: 
        print(f"\n{var} raw values of largest relative {stat} values:")
    else:
        print(f"\n{var} raw values of largest {stat} values:")
    raw_vals = summa.sel(stat=stat0, hru = hru_big)
    # Print all the raw values of the largest nBig values
    for var in plot_vars:
        print(f"{var}: {[f'{val:.1e}' for val in raw_vals[var].values]}")

    # Open the netCDF file with local attributes
    attr = xr.open_dataset(attr_fil)

    # Mask the HRU variable from the netCDF file
    mask = attr['hruId'].isin(hru_big)

    # Get the vegTypeIndex variable from the netCDF file
    vegType_big = attr['vegTypeIndex'][mask]
    lat_big = attr['latitude'][mask]
    lon_big = attr['longitude'][mask]

    # Print the HRU, vegTypeIndex, and values of the largest nBig values
    print("HRU values of the largest values:", hru_big)
    print("vegTypeIndex:", vegType_big.values)
    print("latitude:", np.round(lat_big.values,2))
    print("longitude:", np.round(lon_big.values,2))
    print("Largest error values:", [f"{val:.2e}" for val in val_big.values])