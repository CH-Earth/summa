# written by A. Van Beusekom (2024)

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
    stat = 'mean'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1en'] #cm','be1en','be1lu'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    method_name=['be1','be1en','be1lu'] #maybe make this an argument

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
viz_fl2 = method_name.copy()
for i, m in enumerate(method_name):
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['balance']))
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))

# Specify variables of interest
plot_vars = ['scalarTotalET','scalarCanopyWat','scalarSWE','scalarTotalSoilWat','wallClockTime']
plot_vars2 =['scalarCanopyWat','scalarSWE','scalarTotalSoilWat','averageRoutedRunoff']
comp_vars = ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','numberFluxCalc']
comp_vars2 =['balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass']

plt_titl =  ['(a) Canopy Air Space Energy Balance','(b) Vegetation Energy Balance','(c) Snow Energy Balance','(d) Soil Energy Balance', '(e) Wall Clock Time',]
plt_titl2 = ['(a) Vegetation Mass Balance','(b) Snow Mass Balance','(c) Soil Mass Balance', '(d) Aquifer Mass Balance']
leg_titl0 = ['$W~m^{-3}$'] * 4 +['$num$']
leg_titl02 =['$kg~m^{-2}~s^{-1}$'] * 4

leg_titl = ['mm~y^{-1}$','$kg~m^{-2}$','$kg~m^{-2}$','$kg~m^{-2}$','$s$']
leg_word = ['Total evapotranspiration', 'Total water on the vegetation canopy','Snow water equivalent','Total soil water content','Wall clock time']
leg_titl2 =['$kg~m^{-2}$','$kg~m^{-2}$','$kg~m^{-2}$','mm~y^{-1}$']
leg_word2 = ['Total water on the vegetation canopy','Snow water equivalent','Total soil water content','Average routed runoff']


#fig_fil = '{}_hrly_diff_scat_{}_{}_compressed.png'
#fig_fil = fig_fil.format(','.join(method_name),','.join(settings),stat)
fig_fil = 'BalanceNrg_scat_{}_compressed.png'
fig_fil = fig_fil.format(stat)
fig_fil2 ='BalanceMass_scat_{}_compressed.png'
fig_fil2 =fig_fil2.format(stat)

summa = {}
summa1 = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    summa1[m] = xr.open_dataset(viz_dir/viz_fil[i])
    
##Figure

def run_loop(i,var,comp,leg_t,leg_t0,plt_t,leg_w):
    r = i//2
    c = i-r*2
    if stat == 'mean': 
        word = ' mean'
        stat1 = 'rmnz'
        word1 = ' relative RMSE'
    if stat == 'amax':
        word = ' max'
        stat1 = 'maxe'
        word1 = ' max abs error'
    if comp == 'numberFluxCalc':
        stat1 = stat
        word1 = word

    # Data
    for m in method_name:
        # Get the statistics, remove 9999 (should be nan, but just in case)
        s0 = np.fabs(summa[m][comp].sel(stat=stat)).where(lambda x: x != 9999)
        s = np.fabs(summa1[m][var].sel(stat=stat1)).where(lambda x: x != 9999)

        axs[r,c].scatter(x=s.values,y=s0.values,s=10,zorder=0,label=m)        

    if comp == 'numberFluxCalc':
        stat0_word = 'Number flux calculations'
    else:
        stat0_word = 'Balance abs value'

    stat_word = leg_w
 
    lgnd = axs[r,c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    #axs[r,c].set_title(plt_t)
    #axs[r,c].set_xscale('log')
    if comp != 'numberFluxCalc': axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(stat_word  + word1 + ' [{}]'.format(leg_t))
    axs[r,c].set_ylabel(stat0_word + word + ' [{}]'.format(leg_t0))

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
    
for i,(var,comp,leg_t,leg_t0,plt_t,leg_w) in enumerate(zip(plot_vars,comp_vars,leg_titl,leg_titl0,plt_titl,leg_word)): 
    run_loop(i,var,comp,leg_t,leg_t0,plt_t,leg_w)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)


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
    
for i,(var,comp,leg_t,leg_t0,plt_t,leg_w) in enumerate(zip(plot_vars2,comp_vars2,leg_titl2,leg_titl02,plt_titl2,leg_word2)): 
    run_loop(i,var,comp,leg_t,leg_t0,plt_t,leg_w)

# Save
plt.savefig(viz_dir/fig_fil2, bbox_inches='tight', transparent=False)
