# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python hist_per_GRU.py [stat]
# where stat is rmse or maxe or kgem

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy

viz_dir = Path('/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics')

testing = False
if testing: 
    stat = 'kgem'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be64','be1','sundials_1en6'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    method_name=['be64','be32','be1','sundials_1en6'] #maybe make this an argument

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

#fig_fil = '{}_hrly_diff_hist_{}_{}_zoom_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'Hrly_diff_hist_{}_{}_zoom_compressed.png'
fig_fil = fig_fil.format(','.join(settings),stat)
# possibly want to use these to shrink the axes a bit
if stat=='rmse': maxes = [2,15,8e-6,0.08,7e-9,13e-3]
if stat=='maxe' : maxes = [20,30,3e-4,2,4e-7,0.7]
if stat=='kgem' : maxes = [0.9,0.7,0.9,0.95,0.95,13e-3]

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

    
def run_loop(i,var,mx):
    r = i//2
    c = i-r*2
    num_bins = 200
    stat0 = stat
    if var == 'wallClockTime':
        if stat == 'rmse' or stat == 'kgem': stat0 = 'mean'
        if stat == 'maxe': stat0 = 'amax'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx
    else:
        mx = 0.0
        mn = 1.0
        for m in method_name:
            s = summa[m][var].sel(stat=stat0)
            if stat == 'maxe': s = np.fabs(s) # make absolute value norm
            mx = max(s.max(),mx)
            if stat=='kgem' : mn = min(s.min(),mn)

    # Data
    for m in method_name:
        s = summa[m][var].sel(stat=stat0)
        if stat == 'maxe': s = np.fabs(s) # make absolute value norm
        range = (0,mx)
        if stat=='kgem' and var!='wallClockTime' : range = (mn,1)
        s.plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=range)
    
    if stat == 'rmse': stat_word = ' Hourly RMSE'
    if stat == 'maxe': stat_word = ' Hourly max abs error'
    if stat == 'kgem': stat_word = ' Hourly KGEm'
        
    # wall Clock doesn't do difference
    if var == 'wallClockTime':
        if stat == 'rmse' or stat == 'kgem': stat_word = ' Hourly mean'
        if stat == 'maxe': stat_word = ' Hourly max'

    axs[r,c].legend()
    axs[r,c].set_title(plt_titl[i] + stat_word)
    axs[r,c].set_xlabel('[{}]'.format(leg_titl[i]))
    axs[r,c].set_ylabel('GRU count')
    

for i,(var,mx) in enumerate(zip(plot_vars,maxes)): 
    run_loop(i,var,mx)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
