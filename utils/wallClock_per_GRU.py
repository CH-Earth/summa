# written by A. Van Beusekom (2023)

## Visualize batch number and wallClocktime per GRU, (could be modified to visualize other attributes or forcing per GRU)
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python wallClock_per_GRU.py

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy

viz_dir = Path('/home/avanb/scratch/statistics')

testing = True
if testing: 
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1','be64','sundials_1en6'] #maybe make this an argument
else:
    import sys
    method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))

# Specify variables of interest
plt_titl = ['(a) Wall clock mean time','(b) Wall clock max time']
leg_titl = ['$s$','$s$']

#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'WallClockTime_batchNum_scat_compressed.png'

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
    fig,axs = plt.subplots(1,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(1,2,figsize=(140,133))

    
def run_loop(i,stat):
    r = i//2
    c = i-r*2

    # Data
    for m in method_name:
        s = summa[m]['wallClockTime'].sel(stat=stat)
        modulus = s.indexes['hru'] % 518
        axs[c].scatter(x=s.values,y=modulus.values,s=1,zorder=0,label=m)
        
    if stat == 'mean': stat_word = 'Wall clock time hourly mean '
    if stat == 'amax': stat_word = 'Wall clock time hourly max '
    stat0_word ='Batch number'

 
    lgnd = axs[c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[c].set_title(plt_titl[i])
    axs[c].set_xlabel(stat_word + '[{}]'.format(leg_titl[i]))
    axs[c].set_ylabel(stat0_word)


for i,stat in enumerate(['mean','amax']): 
    run_loop(i,stat)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)

