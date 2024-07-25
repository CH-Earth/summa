# written by A. Van Beusekom (2024)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python wallStd_per_GRU.py

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy
import pandas as pd

run_local = False
if run_local: 
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1en','be1en'] #cm','be1en','be1lu'] #maybe make this an argument
else:
    import sys
    viz_dir = Path('/home/avanb/scratch/statistics')
    method_name=['sundials_1en4','sundials_1en6','sundials_1en8'] #maybe make this an argument

# Simulation statistics file locations
viz_fl2 = method_name.copy()
for i, m in enumerate(method_name):
    viz_fl2[i] = m + '_hrly_diff_wals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['wallclock']))

# Specify variables of interest
stats = ['mean','amax']
plot_vars = ['wallClockTime']*2
comp_vars = ['wallClockTime']*2

plt_titl =  ['(a) Wall Clock Time Mean','(b) Wall Clock Time Max']
leg_titl0 = ['$s$']*2
leg_titl = ['$s$']*2

#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'WallStd_scat_mean_max_compressed.png'

summa = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    
##Figure

def run_loop(i,var,comp,leg_t,leg_t0,plt_t,stat):
    r = i//2
    c = i-r*2

    # Data
    for mm,m in enumerate(method_name):
        # Get the statistics, remove 9999 (should be nan, but just in case)
        s = np.fabs(summa[m][var].sel(stat=stat)).where(lambda x: x != 9999)
        s0 = np.fabs(summa[m][comp].sel(stat='std')).where(lambda x: x != 9999)

        stat_word = 'Time'
        stat0_word = 'Time std dev'

        axs[c].scatter(x=s.values,y=s0.values,s=10,zorder=0,label=m)

        # Create a mask that is True where `s.values` and `s0.values` are not NaN
        mask = ~np.isnan(s.values) & ~np.isnan(s0.values)

        # Use the mask to filter `s.values` and `s0.values`
        filtered_s_values = s.values[mask]
        filtered_s0_values = s0.values[mask]

        # Fit a linear regression model
        coefficients = np.polyfit(filtered_s_values, filtered_s0_values, 1)
        polynomial = np.poly1d(coefficients)

        # Calculate the R-squared value
        correlation_matrix = np.corrcoef(filtered_s_values, filtered_s0_values)
        correlation_xy = correlation_matrix[0,1]
        r_squared = correlation_xy**2
        # Add the R-squared value to the plot
        axs[c].text(0.8, 0.5-0.03*mm, f'{m} R² = {r_squared:.2f}', transform=axs[c].transAxes)
        print(m,stat,'Coefficients:', coefficients, 'R-squared:', r_squared)


    if stat == 'mean': word = ' mean'
    if stat == 'amax': word = ' max'
 
    lgnd = axs[c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[c].set_title(plt_t)
    axs[c].set_xscale('log')
    axs[c].set_yscale('log')
    
    axs[c].set_xlabel(stat_word  + word + ' [{}]'.format(leg_t))
    axs[c].set_ylabel(stat0_word + ' [{}]'.format(leg_t0))

if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(2,figsize=(35,38))
else:
    fig,axs = plt.subplots(2,figsize=(140,160))
fig.subplots_adjust(hspace=0.24, wspace=0.24) # Adjust the bottom margin, vertical space, and horizontal space
#fig.suptitle('Scatterplot of Hourly Statistics for each GRU', fontsize=40,y=1.0)
    
for i,(var,comp,leg_t,leg_t0,plt_t,stat) in enumerate(zip(plot_vars,comp_vars,leg_titl,leg_titl0,plt_titl,stats)): 
    run_loop(i,var,comp,leg_t,leg_t0,plt_t,stat)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)

