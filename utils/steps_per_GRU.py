# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python steps_per_GRU.py [stat]
# where stat is mean or amax

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
    stat = 'amax'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')
    method_name=['be1','be4','be8','be16','be32'] #maybe make this an argument
    method_name=['be1'] #maybe make this an argument
    #method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
stepsets= ['numberStateSplit','numberDomainSplitNrg','numberDomainSplitMass','numberScalarSolutions','meanStepSize']

viz_fil = method_name.copy()
viz_fl2 = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
    viz_fl2[i] = m + '_hrly_diff_steps_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(stepsets))

# Specify variables of interest
plot_vars = ['numberDomainSplitNrg','numberDomainSplitMass','numberScalarSolutions','meanStepSize']
plt_titl = ['(a) Energy Domain Splits','(b) Mass Domain Splits','(c) Scalar Solutions', '(d) Mean Step Size']
leg_titl = ['$num$', '$num$','$num$','$s$']
leg_titl0 = ['$num$', '$num$','$num$','$s$']

#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'Splits_steps_scat_{}_compressed.png'
fig_fil = fig_fil.format(stat)

summa = {}
wall = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    wall[m] = xr.open_dataset(viz_dir/viz_fil[i])
    
##Figure

# Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(2,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(2,2,figsize=(140,133))
fig.suptitle('Hourly Splits and Time Steps  for each GRU', fontsize=40)
    
def run_loop(i,var):
    r = i//2
    c = i-r*2

    # Data
    for m in method_name:
        s = summa[m][var].sel(stat=stat)

        if var == 'numberDomainSplitNrg' or var == 'numberDomainSplitMass':
            s0 = summa[m]['numberStateSplit'].sel(stat=stat)
            stat0_word = 'State Splits'
            if var == 'numberDomainSplitNrg': stat_word = 'Energy Domain Splits'
            if var == 'numberDomainSplitMass': stat_word = 'Mass Domain Splits'
        if var == 'numberScalarSolutions':
            s0 = summa[m]['numberDomainSplitMass'].sel(stat=stat) + summa[m]['numberDomainSplitNrg'].sel(stat=stat)
            stat0_word = 'Domain Splits'
            stat_word = 'Scalar Solutions'
        if var == 'meanStepSize':
            s0 = wall[m]['wallClockTime'].sel(stat=stat)
            stat0_word = 'Wallclock Time'
            stat_word = 'Mean Step Size'

        axs[r,c].scatter(x=s.values,y=s0.values,s=10,zorder=0,label=m)        
 
    if stat == 'mean': word = ' mean'
    if stat == 'amax': word = ' max'
 
    lgnd = axs[r,c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[r,c].set_title(plt_titl[i])
    axs[r,c].set_xlabel(stat_word  + word + ' [{}]'.format(leg_titl[i]))
    axs[r,c].set_ylabel(stat0_word + word + ' [{}]'.format(leg_titl0[i]))


for i,var in enumerate(plot_vars): 
    run_loop(i,var)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
