# written originally by W. Knoben, modified by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    Catchment shapefile with GRU delineation
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# To improve visualization of large lakes, HydroLAKES lake delineations are plotted on top of the catchment GRUs and river segments.
# Dealing with HydroLAKES inputs is not considered within scope of the workflow and therefore requires some manual downloading and preprocessing of this data for those who wish to reproduce this step.
# The relevant code is easily disabled by switching the plot_lakes = True flag to False.

# Run:
# python hist_per_GRU.py rmse


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
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    stat = "max" 
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1] # max or rmse
 
# Simulation statistics file locations
method_name=['sundials_1en6', 'be1','be32'] #maybe make this an argument

settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}
viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings.keys()),','.join(settings.values()))

# Specify variables of interest
plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}~s^{-1}$','$kg~m^{-2}$','$m~s^{-1}$','$s$']

fig_fil =''
for m in method_name:
    fig_fil = fig_fil + m
fig_fil = fig_fil + '_hrly_diff_hist_{}_{}__zoom_compressed.png'
fig_fil = fig_fil.format(','.join(settings.keys()),','.join(settings.values()))
# possibly want to use these to shrink the axes a bit
if stat=='rmse': maxes = [2,15,8e-6,0.08,7e-9,7e-3]
if stat=='max' : maxes = [20,30,3e-4,2,35e-8,0.7]


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
    if 'zoom' in fig_fil:
        mx = mx
    else:
        mx = 0.0
        for m in method_name:
            s = summa[m]
            mx = max(s[var].max(),mx)

    # Data
    for m in method_name:
        s = summa[m]
        s[var].plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=(0,mx))
    
    if stat == 'rmse':
        stat_word = ' Hourly RMSE'
    if stat == 'max':
        stat_word = ' Hourly max abs error'
        
    # wall Clock doesn't do difference
    if var == 'wallClockTime':
        if stat == 'rmse': stat_word = ' Hourly mean'
        if stat == 'max': stat_word = ' Hourly max'

    axs[r,c].legend()
    axs[r,c].set_title(plt_titl[i] + stat_word)
    axs[r,c].set_xlabel('[{}]'.format(leg_titl[i]))
    axs[r,c].set_ylabel('GRU count')
    

for i,(var,mx) in enumerate(zip(plot_vars,maxes)): 
    run_loop(i,var,mx)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)
