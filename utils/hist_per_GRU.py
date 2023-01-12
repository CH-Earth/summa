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
method_name1='sundials_1en6' #keep constant?
method_name2='be1' #maybe make this an argument

settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}
viz_fil1 = method_name1 + '_hrly_diff_stats_{}_{}.nc'
viz_fil1 = viz_fil1.format(','.join(settings.keys()),','.join(settings.values()))
viz_fil2 = method_name2 + '_hrly_diff_stats_{}_{}.nc'
viz_fil2 = viz_fil2.format(','.join(settings.keys()),','.join(settings.values()))


# Specify variables of interest
plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}~s^{-1}$','$kg~m^{-2}$','$m~s^{-1}$','$s$']

fig_fil = method_name1 + method_name2 + '_hrly_diff_hist_{}_{}_compressed.png'
fig_fil = fig_fil.format(','.join(settings.keys()),','.join(settings.values()))
# possibly want to use these to shrink the axes a bit
if stat=='rmse': maxes = [120,80,3e-5,0.2,1e-8,1e-2]
if stat=='max' : maxes = [100,50,5e-4,4,35e-8,100]

## Pre-processing, map SUMMA sims to catchment shapes
# Get the aggregated statistics of SUMMA simulations
summa1 = xr.open_dataset(viz_dir/viz_fil1)
summa2 = xr.open_dataset(viz_dir/viz_fil2)

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
    vmin = 0
    num_bins = 200
    mx = max([summa1[var].max(),summa2[var].max()]).values

    # Data
    sm = summa1[var].plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=method_name1,linewidth=2.0,range=(0,mx))
    sm = summa2[var].plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=1,label=method_name2,linewidth=2.0,range=(0,mx))
    
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
