#!/usr/bin/env python

# python largest_error_attrib.py [method_name] [stat]

import numpy as np
import xarray as xr
from pathlib import Path

nBig = 10
do_rel = False # stat relative to the benchmark simulation
do_var = True # do vars, if False do bals

run_local = True
if run_local:
    top_fold = '/Users/amedin/Research/USask/test_py/'
    attr_fold = '/Users/amedin/Research/USask/test_py/settings/'
    method_name= 'be8en'
    stat = 'avge'
else:
    import sys
    top_fold    = '/home/avanb/scratch/'
    attr_fold = '/home/avanb/TestScripts/settings/'
    method_name = sys.argv[1]
    stat = sys.argv[2]

des_dir =  top_fold + 'statistics_en'
des_dir = Path(des_dir)

if do_var:
    settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','scalarRootZoneTemp']
    settings= ['scalarTotalSoilWat','scalarRootZoneTemp']
    viz_fil = method_name + '_hrly_diff_stats_{}.nc'
    viz_fil = viz_fil.format(','.join(['accuracy']))
    src_file =  des_dir / viz_fil
    plot_vars = settings.copy()
    short_name= [#'SWE     ',
                 'soilWat ',
                 #'ET      ',
                 #'canWat  ',
                 'rootTemp']
else:
    do_rel = False
    settings= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']
    viz_fil = method_name + '_hrly_diff_bals_balance.nc'
    src_file =  des_dir / viz_fil
    plot_vars = settings.copy()
    short_name= ['casNrg  ',
                 'vegNrg  ',
                 'snowNrg ',
                 'soilNrg ',
                 'vegMass ',
                 'snowMass',
                 'soilMass',
                 'aqMass  ']

attr_fil = Path(attr_fold) / 'attributes.nc'

# Open the netCDF file with RMSE data
summa = xr.open_dataset(src_file)
if stat == 'rmse' or stat == 'kgem' or stat == 'mean' or stat == 'avge': statr = 'mean_ben'
if stat == 'rmnz' or stat == 'mnnz': statr = 'mnnz_ben'
if stat == 'maxe' or stat == 'amax': statr = 'amax_ben'

for var in plot_vars:

    # Get the variable from the netCDF file
    stat0 = stat
    if var == 'wallClockTime': 
        if stat == 'rmse' or stat == 'kgem' or stat == 'mean' or stat == 'avge': stat0 = 'mean'
        if stat == 'rmnz' or stat == 'mnnz': stat0 = 'mnnz'
        if stat == 'maxe' or stat == 'amax': stat0 = 'amax'
    s = summa[var].sel(stat=stat0)
    if do_rel: 
        s_rel = summa[var].sel(stat=statr)
        if var != 'wallClockTime': s = s/s_rel

    # Mask the finite values of the variable
    mask = np.isfinite(s)
    s = s[mask]

    # Get the indices of the largest nBig values
    big_indices = abs(s).argsort()[-nBig:]
    
    # Get the largest nBig values
    val_big = s[big_indices.values]
    
    # Get the hru coordinate of the largest nBig values
    hru_big = s[big_indices.values].hru.values #gru index

    # Get the largest nBig bench values
    if do_rel: ben_big = s_rel.sel(hru=hru_big)
    
    # Print all the values of the biggest rmse hru
    if do_rel: 
        print(f"\n{var} raw error values of largest relative {stat} values:")
    else:
        print(f"\n{var} raw error values of largest {stat} values:")
    # Print all the raw values of the largest nBig values
    raw_vals = summa.sel(stat=stat0, hru=hru_big)
    if do_var: 
        plot_var = plot_vars
    else: 
        plot_var = plot_vars[:-1]
    for i,var0 in enumerate(plot_var):
        print(f"{short_name[i]}: [{' '.join(f'{val:8.1e}' for val in raw_vals[var0].values)}]")
    var0 = 'wallClockTime'
    if stat == 'rmse' or stat == 'kgem' or stat == 'mean': stat00 = 'mean'
    if stat == 'rmnz' or stat == 'mnnz': stat00 = 'mnnz'
    if stat == 'maxe' or stat == 'amax': stat00 = 'amax'
    raw_vals = summa.sel(stat=stat00, hru=hru_big)
    if not do_var: print("wall"f"{stat00}: [{' '.join(f'{val:8.1e}' for val in raw_vals[var0].values)}]")

    # Open the netCDF file with local attributes
    attr = xr.open_dataset(attr_fil)

    # Mask the HRU variable from the netCDF file
    #mask = attr['hruId'].isin(hru_big)

    # Filtered HRU IDs
    #filtered_hru_ids = attr['hruId'][mask]
    
    # Determine the indices that would sort filtered_hru_ids to match the order of hru_big
    #h_ind = [filtered_hru_ids.values.tolist().index(hru_id) for hru_id in hru_big if hru_id in filtered_hru_ids.values]
    #h = attr['hruId'][mask].values[h_ind]

    # Get the vegTypeIndex, lat, lon variables from the netCDF file
    #vegType_big = attr['vegTypeIndex'][mask].values[h_ind]
    #lat_big = attr['latitude'][mask].values[h_ind]
    #lon_big = attr['longitude'][mask].values[h_ind]


    # with Actors
    filtered_hru_ids = attr['hruId'][hru_big-1].values

    # Get the vegTypeIndex, lat, lon variables from the netCDF file
    vegType_big = attr['vegTypeIndex'][hru_big-1].values
    lat_big = attr['latitude'][hru_big-1].values
    lon_big = attr['longitude'][hru_big-1].values

    # Print the attributes of the largest nBig values
    print("HRU vals: [", " ".join([f"{val:8d}"  for val in filtered_hru_ids]), "]", sep="")
    print("vegType : [", " ".join([f"{val:8d}"  for val in vegType_big]), "]", sep="")
    print("latitude: [", " ".join([f"{val:8.2f}"  for val in lat_big]), "]", sep="")
    print("longitud: [", " ".join([f"{val:8.2f}"  for val in lon_big]), "]", sep="")

    # Print the values of the largest nBig values, bench will be the mean, mnnz, or amax and err will be the rmse, rmnz, or maxe
    if do_rel: 
        #print(summa[var].sel(stat=stat0, hru=hru_big))
        print("Ben vals: [", " ".join([f"{val:8.1e}" for val in ben_big.values]), "]", sep="")
    print("Err vals: [", " ".join([f"{val:8.1e}" for val in val_big.values]), "]", sep="")
    