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
# python hist_per_GRU.py sundials_1en6 rmse


# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy

#viz_dir = Path('/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics')
viz_dir = Path('/home/avanb/scratch/statistics')

testing = False
if testing: 
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name = 'be1'
    stat = "rmse" 
else:
    import sys
    # The first input argument specifies the run where the files are
    method_name = sys.argv[1] # sys.argv values are strings by default so this is fine (sundials_1en6 or be1)
    stat = sys.argv[2] # max or rmse
 
# Simulation statistics file locations
settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}
viz_fil = method_name + '_hrly_diff_stats_{}_{}.nc'
viz_fil = viz_fil.format(','.join(settings.keys()),','.join(settings.values()))


# Specify variables of interest
plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}~s^{-1}$','$kg~m^{-2}$','$m~s^{-1}$','$s$']

fig_fil = method_name + '_hrly_diff_hist_{}_{}_compressed.png'
fig_fil = fig_fil.format(','.join(settings.keys()),','.join(settings.values()))

## Pre-processing, map SUMMA sims to catchment shapes
# Get the aggregated statistics of SUMMA simulations
summa = xr.open_dataset(viz_dir/viz_fil)

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

plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

# colorbar axes
f_x_mat = [0.473,0.97,0.473,0.97,0.473,0.97]
f_y_mat = [0.69,0.69,0.36,0.36,0.03,0.03]

#plt.tight_layout()

    
def run_loop(i,var,f_x,f_y):
    r = i//2
    c = i-r*2
    
    # Data
    sm = summa[var].plot.hist(ax=axs[r,c])
    
    # Custom colorbar
    #cax = fig.add_axes([f_x,f_y,0.02,0.3])
    #cbr = fig.colorbar(sm, cax=cax, extend='max')
    #cbr.ax.set_title('[{}]'.format(leg_titl[i]))
    
    if stat == 'rmse':
        stat_word = ' Hourly RMSE'
    if stat == 'max':
        stat_word = ' Hourly max abs error'
        
    # wall Clock doesn't do difference
    if var == 'wallClockTime':
        if stat == 'rmse': stat_word = ' Hourly mean'
        if stat == 'max': stat_word = ' Hourly max'

    axs[r,c].set_title(plt_titl[i] + stat_word)
    axs[r,c].set_xlabel('[{}]'.format(leg_titl[i]))
    

for i,(var,f_x,f_y) in enumerate(zip(plot_vars,f_x_mat,f_y_mat)): 
    run_loop(i,var,f_x,f_y)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)
