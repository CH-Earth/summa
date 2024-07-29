# Combine split domain state files (with 2 dimensions, hru and gru)
# Modified by W. Knoben (2021) from A. Wood (2020)
# -----------------------------------------------

import sys, glob
import pandas as pd
import xarray as xr

# --------- arguments -----------
'''
print("found %d args" % len(sys.argv))
if len(sys.argv) == 5:
    stateFileRoot = sys.argv[1]     # eg './hstate/wbout_restart_'
    startDate     = sys.argv[2]     # eg '19910301'
    endDate       = sys.argv[3]     # eg '19920301'
    Freq          = sys.argv[4]     # D (daily) or MS (monthly)
else:
    print("USAGE: %s input_filepath/root startdate(YYYYMMDD) enddate frequency_of_states(D or MS)" % sys.argv[0])
    sys.exit(0)
'''
srcPath = '/project/gwf/gwf_cmt/wknoben/summaWorkflow_data/domain_Nelson/simulations/run3_be4_make_ics/SUMMA'
srcName = 'run3_be4_make_ics_restart_2017123123_*.nc'
desPath = '/project/gwf/gwf_cmt/wknoben/summaWorkflow_data/domain_Nelson/settings/SUMMA/'
desName = 'warmState.nc'

# --------- code -----------
# find the files
output_file_list = glob.glob(srcPath + '/' + srcName)
output_file_list.sort()

out_ds   = [xr.open_dataset(f) for f in output_file_list]
hru_vars = [] # variables that have hru dimension
gru_vars = [] # variables that have gru dimension

for name, var in out_ds[0].variables.items():
    if 'hru' in var.dims:
        hru_vars.append(name)
    elif 'gru' in var.dims:
        gru_vars.append(name)

hru_ds = [ds[hru_vars] for ds in out_ds]
gru_ds = [ds[gru_vars] for ds in out_ds]
hru_merged = xr.concat(hru_ds, dim='hru')
gru_merged = xr.concat(gru_ds, dim='gru')
merged_ds  = xr.merge([hru_merged, gru_merged])

merged_ds.load().to_netcdf(desPath + '/' + desName)