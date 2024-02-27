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

viz_dir = Path('/home/avanb/scratch/statistics')

testing = False
if testing: 
    stat = 'amax'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    method_name=['be1','be1en','be1ln'] #maybe make this an argument

# Simulation statistics file locations
balssets= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']

viz_fil = method_name.copy()
viz_fl2 = method_name.copy()
for i, m in enumerate(method_name):
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(balssets))

# Specify variables of interest
plot_vars = ['balanceVegNrg','balanceSnowNrg','balanceSoilNrg','wallClockTime','balanceCasNrg']
plot_vars2 = ['balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass']

plt_titl = ['(a) Vegetation Energy Balance','(b) Snow Energy Balance','(c) Soil Energy Balance', '(d) Wall Clock Time','e) Canopy Air Space Energy Balance']
plt_titl2 = ['(a) Vegetation Mass Balance','(b) Snow Mass Balance','(c) Soil Mass Balance', '(d) Aquifer Mass Balance']
leg_titl = ['$W~m^{-2}$', '$W~m^{-2}$','$W~m^{-2}$','$s$','$W~m^{-2}$']
leg_titl2 = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}$','$kg~m^{-2}$']
leg_titl0 = ['$W~m^{-2}$', '$W~m^{-2}$','$W~m^{-2}$','$num$','$W~m^{-2}$']
leg_titl20 = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}$','$kg~m^{-2}$']

#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'Balance_scat_{}_compressed.png'
fig_fil = fig_fil.format(stat)

summa = {}
wall = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    
##Figure

if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(3,2,figsize=(35,38))
else:
    fig,axs = plt.subplots(3,2,figsize=(140,160))
fig.subplots_adjust(hspace=0.24, wspace=0.24) # Adjust the bottom margin, vertical space, and horizontal space
#fig.suptitle('Scatterplot of Hourly Statistics for each GRU', fontsize=40,y=1.0)
    
def run_loop(i,var):
    r = i//2
    c = i-r*2

    # Data
    for m in method_name:
        s = summa[m][var].sel(stat='mean')

        if var == 'wallClockTime':
            s0 = summa[m]['numberFluxCalc'].sel(stat='mean')
            stat0_word = 'Mean Number Flux Calculations'
            stat_word = 'Mean Wallclock Time'
        else:
            s0 = summa[m]['var'].sel(stat='amax')
            stat0_word = 'Max Absolute Value'
            stat_word = 'Mean Absolute Value'

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
