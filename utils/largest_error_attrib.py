#!/usr/bin/env python

# python largest_error_attrib.py [method_name] 

import sys
import numpy as np
import xarray as xr
from pathlib import Path

top_fold    = '/home/avanb/scratch/'
attr_fold = '/home/avanb/TestScripts/settings/'
nBig = 10

testing = False
if testing:
    method_name= 'be1' #sys.argv[1]
    top_fold = '/Users/amedin/Research/USask/test_py/'
    attr_fold = '/Users/amedin/Research/USask/test_py/settings/'
    nBig = 10
else:
    method_name = sys.argv[1]

des_dir =  top_fold + 'statistics'
des_dir = Path(des_dir)

settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name + '_hrly_diff_stats_{}.nc'
viz_fil = viz_fil.format(','.join(settings))
src_file =  des_dir / viz_fil

attr_fil = Path(attr_fold) / 'attributes.nc'

# Open the netCDF file with RMSE data
summa = xr.open_dataset(src_file)

# Get the RMSE variable from the netCDF file
rmse_var = summa['scalarTotalSoilWat'].sel(stat='rmse')
#rmse_var = summa['scalarTotalET'].sel(stat='rmse')
#rmse_var = summa['averageRoutedRunoff'].sel(stat='rmse')

# Mask the finite values of the RMSE variable
mask = np.isfinite(rmse_var)
rmse_var = rmse_var[mask]

# Get the indices of the largest nBig RMSE values
max_rmse_indices = rmse_var.argsort()[-nBig:]

# Get the largest RMSE values
rmse_of_max_rmse = rmse_var[max_rmse_indices.values]

# Get the hru coordinate of the largest RMSE values
hru_of_max_rmse = rmse_var[max_rmse_indices.values].hru.values

# Print all the values of the the biggest rmse hru
rmse_values = summa.sel(hru=hru_of_max_rmse[nBig-1],stat='rmse')
print(settings[0:5])
print(rmse_values[settings[0]].values,rmse_values[settings[1]].values,rmse_values[settings[2]].values,rmse_values[settings[3]].values,rmse_values[settings[4]].values)

# Open the netCDF file with local attributes
attr = xr.open_dataset(attr_fil)

# Mask the HRU variable from the netCDF file
mask = attr['hruId'].isin(hru_of_max_rmse)

# Get the vegTypeIndex variable from the netCDF file
veg_type_of_max_rmse = attr['vegTypeIndex'][mask]

# Print the HRU, vegTypeIndex, and RMSE values of the largest RMSE values
print("HRU values of the largest RMSE values:", hru_of_max_rmse)
print("GRU values of the largest RMSE values:", gru_of_max_rmse)
print("vegTypeIndex values of the largest RMSE values:", veg_type_of_max_rmse.values)
print("RMSE values of the largest RMSE values:", rmse_of_max_rmse.values)