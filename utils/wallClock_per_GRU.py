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
import pandas as pd

run_local = False
if run_local: 
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1'] #maybe make this an argument
else:
    import sys
    viz_dir = Path('/home/avanb/scratch/statistics')
    method_name=['be1','be4','be8','be16','be32'] #sundials will not show node differences as much

nbatch_hrus = 518 # number of HRUs per batch
# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
eff_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
    eff_fil[i] = 'eff_' + m + '.txt'

# Specify variables of interest
plt_titl = ['(a) Mean time vs basin in batch','(b) Max time vs basin in batch','(c) Mean time vs node','(d) Max time vs node']
leg_titl = ['$s$','$s$','$s$','$s$']

#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'WallClockTime_batchNum_scat_compressed.png'

summa = {}
eff = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
    # Read the data from the eff.txt file into a DataFrame
    eff[m] = pd.read_csv(viz_dir/eff_fil[i], sep=',', header=None, names=['CPU Efficiency', 'Array ID', 'Job Wall-clock time', 'Node Number'])
    # Extract only the values after the ':' character in the 'CPU Efficiency', 'Job Wall-clock time', and 'Node Number' columns
    eff[m]['CPU Efficiency'] = eff[m]['CPU Efficiency'].str.split(':').str[1].astype(float)
    eff[m]['Array ID'] = eff[m]['Array ID'].str.split(':').str[1].astype(int)   
    eff[m]['Job Wall-clock time'] = eff[m]['Job Wall-clock time'].str.split(':').str[1].astype(float)
    eff[m]['Node Number'] = eff[m]['Node Number'].str.split(':').str[1].astype(int)

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
    fig,axs = plt.subplots(2,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(2,2,figsize=(140,133))
fig.suptitle('BE Hourly Wallclock for each GRU', fontsize=40)
    
def run_loop(c,stat):

    # Data
    for m in method_name:
        s = summa[m]['wallClockTime'].sel(stat=stat)
        batch = np.floor(np.arange(len(s.indexes['hru'])) /nbatch_hrus)
        basin_num = np.arange(len(s.indexes['hru'])) % nbatch_hrus
        # Create a dictionary to store the values for each batch
        efficiency = {}
        node = {}
        # Iterate over the rows in the data DataFrame
        for index, row in eff[m].iterrows():
            # Extract the values from the row
            batch0 = int(row['Array ID'])
            eff0 = row['CPU Efficiency']
            node0 = row['Node Number']
            # Store the value for the current batch in the dictionary
            efficiency[batch0] = eff0  
            node[batch0] = node0
        # Select the values for the current batch using boolean indexing
        eff_batch = np.array([efficiency[b] for b in batch])
        node_batch = np.array([node[b] for b in batch])
        # Multiply the efficiency values by the s values
        x = eff_batch * s.values

        if stat == 'mean': stat_word = 'Wall clock time hourly mean * efficiency '
        if stat == 'amax': stat_word = 'Wall clock time hourly max * efficiency '
        for r in range(2):
            if r == 0:
                axs[r,c].scatter(x=x,y=basin_num,s=1,zorder=0,label=m)
                stat0_word ='Basin number in batch'
                lgnd0 = axs[r,c].legend()
            if r == 1:
                axs[r,c].scatter(x=x,y=node_batch,s=1,zorder=0,label=m)
                stat0_word ='Node number'
                lgnd1 = axs[r,c].legend()
            axs[r,c].set_title(plt_titl[r*2+c])
            axs[r,c].set_xlabel(stat_word + '[{}]'.format(leg_titl[r*2+c]))
            axs[r,c].set_ylabel(stat0_word)
    
    for j, m in enumerate(method_name):
        lgnd0.legendHandles[j]._sizes = [80]
        lgnd1.legendHandles[j]._sizes = [80]


for i,stat in enumerate(['mean','amax']): 
    run_loop(i,stat)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
