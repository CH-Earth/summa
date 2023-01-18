# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python scat_per_GRU.py rmse
# where stat is rmse or maxe

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import copy

viz_dir = Path('/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics')

testing = False
if testing: 
    stat = 'rmse'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['sundials_1en6','be1'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1] # max or rmse
    method_name=['sundials_1en6','be32','be1'] #maybe make this an argument

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))

# Specify variables of interest
plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}~s^{-1}$','$kg~m^{-2}$','$m~s^{-1}$','$s$']

fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
# possibly want to use these to shrink the axes a bit
if stat=='rmse': 
    maxes = [2,15,8e-6,0.08,7e-9,100]
    maxes0 = [2,15,8e-6,0.08,7e-9,13e-3]
if stat=='maxe' : 
    maxes = [20,30,3e-4,2,4e-7,5000]
    maxes0 = [100,2000,-1e-3,20,2e-5,0.7]

# Get the aggregated statistics of SUMMA simulations
summa = {}
for i, m in enumerate(method_name):
    summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
    
##Figure

# Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

# Flip the evaporation values so that they become positive, not if plotting diffs
#bas_albers['plot_ET'] = bas_albers['scalarTotalET'] * -1
#bas_albers['plot_ET'] = bas_albers['plot_ET'].where(bas_albers['scalarTotalET'] != -9999, np.nan)

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(3,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(3,2,figsize=(140,133))

    
def run_loop(i,var,mx,mx0):
    r = i//2
    c = i-r*2
    if stat == 'rmse': stat0 = 'mean'
    if stat == 'maxe': stat0 = 'amax'

    # Data
    for m in method_name:
        s = summa[m][var].sel(stat=[stat,stat0])
        if stat == 'maxe': s.loc[dict(stat='maxe')] = np.fabs(s.loc[dict(stat='maxe')]) # make absolute value norm
        axs[r,c].scatter(x=s.sel(stat=stat).values,y=s.sel(stat=stat0).values,zorder=0,label=m)
        if 'zoom' in fig_fil: 
            axs[r,c].set_xlim(0,mx)
            if mx0<0: axs[r,c].set_ylim(mx0,0)
            if mx0>0: axs[r,c].set_ylim(0,mx0)
        
    if stat == 'rmse': 
        stat_word = 'Hourly RMSE '
        stat0_word ='Hourly mean '
    if stat == 'maxe': 
        stat_word = ' Hourly max abs error '
        stat0_word =' Hourly max '
 
    axs[r,c].legend()
    axs[r,c].set_title(plt_titl[i])
    axs[r,c].set_xlabel(stat_word + '[{}]'.format(leg_titl[i]))
    axs[r,c].set_ylabel(stat0_word + '[{}]'.format(leg_titl[i]))
    

for i,(var,mx,mx0) in enumerate(zip(plot_vars,maxes,maxes0)): 
    run_loop(i,var,mx,mx0)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)

