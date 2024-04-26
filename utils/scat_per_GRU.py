# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python scat_per_GRU.py [stat]
# where stat is rmse or maxe or kgem

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
nbatch_hrus = 518 # number of HRUs per batch
do_rel = False # plot relative to the benchmark simulation

# which statistics to plot
do_vars = True
do_balance = True

testing = False
if testing: 
    stat = 'rmnz'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1en']
    plt_name=['BE1 mixed']
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    #method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument
    #plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6'] #maybe make this an argument
    #method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument
    #plt_name=['BE1','BE16','BE32','IDAe-6'] #maybe make this an argument
    method_name=['be1','be1en']
    plt_name=['BE1 temp','BE1 mixed']

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
viz_fl2 = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['balance']))

summa = {}
summa1 = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
    summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    
def run_loop(i,var,plt_t,leg_t,leg_t0,leg_tm):

    r = i//2
    c = i-r*2
    if stat == 'rmse' or stat == 'kgem': 
        stat0 = 'mean'
        stat0_word = 'mean'
        statr = 'mean_ben'
        stat_word = 'RMSE' 
        if stat == 'kgem': stat_word = 'KGE"'

    if stat == 'rmnz':
        stat0 = 'mnnz'
        stat0_word = 'mean' # no 0s'
        statr = 'mnnz_ben'
        stat_word = 'RMSE' # no 0s'
    if stat == 'maxe': 
        stat0 = 'amax'
        stat0_word = 'max'
        statr = 'amax_ben'
        stat_word = 'max abs error'

    # Data
    if do_rel: s_rel = summa[method_name[0]][var].sel(stat=statr)
    for m in method_name:
        s = summa[m][var].sel(stat=[stat,stat0])
        if do_rel and var != 'wallClockTime': s = s/s_rel

        if var == 'scalarTotalET' and not do_rel:
            if stat =='rmse' or stat =='rmnz' : s = s*31557600 # make annual total
            if stat =='maxe': s = s*3600 # make hourly max
        if var == 'averageRoutedRunoff'and not do_rel:
            if stat =='rmse' or stat =='rmnz' : s = s*31557600*1000 # make annual total
            if stat =='maxe': s = s*3600*1000 # make hourly max      
        if stat == 'maxe': s.loc[dict(stat='maxe')] = np.fabs(s.loc[dict(stat='maxe')]) # make absolute value norm

        axs[r,c].scatter(x=np.fabs(s.sel(stat=stat).values),y=s.sel(stat=stat0).values,s=1,zorder=0,label=m)        
 
    lgnd = axs[r,c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[r,c].set_title(plt_t)
    if stat == 'rmse' or stat == 'rmnz': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_t))
    if stat == 'maxe': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_tm))   
    if stat == 'kgem': axs[r,c].set_xlabel(stat_word)
    #if do_rel and var!='wallClockTime': axs[r,c].set_xlabel(stat_word + ' rel to bench ' + stat0_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    axs[r,c].set_ylabel(stat0_word + ' [{}]'.format(leg_t0))
    #if do_rel and var!='wallClockTime': axs[r,c].set_ylabel(stat0_word + ' rel to bench ' + stat0_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_ylabel('relative '+ stat0_word)


def run_loopb(i,var,comp,leg_t,leg_t0,plt_t):
    r = i//2
    c = i-r*2

    if stat == 'rmse' or stat == 'kgem': 
        stat0 = 'mean'
        word = ' mean'
    if stat == 'rmnz':
        stat0 = 'mean'
        word = ' mean' # no 0s'
    if stat == 'maxe': 
        stat0 = 'amax'
        word = ' max'

    # Data
    for m in method_name:
        # Get the statistics, remove 9999 (should be nan, but just in case)
        s0 = np.fabs(summa1[m][comp].sel(stat=stat0)).where(lambda x: x != 9999)
        s = np.fabs(summa1[m][var].sel(stat=stat0)).where(lambda x: x != 9999)

        axs[r,c].scatter(x=s.values,y=s0.values,s=10,zorder=0,label=m)        

    if comp == 'numberFluxCalc':
        stat0_word = 'Number flux calculations'
        stat_word = 'Wall clock time'
    else:
        stat0_word = 'Balance abs value'
        stat_word = 'Balance abs value'
 
    lgnd = axs[r,c].legend()
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[r,c].set_title(plt_t)
    axs[r,c].set_xscale('log')
    if comp != 'numberFluxCalc': axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(stat_word  + word + ' [{}]'.format(leg_t))
    axs[r,c].set_ylabel(stat0_word + word + ' [{}]'.format(leg_t0))


if do_vars:

    plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff']
    leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$']
    leg_titl0 = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$']
    leg_titlm= ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~h^{-1}$','$kg~m^{-2}$','$mm~h^{-1}$']

    fig_fil = 'Hrly_diff_scat_{}_{}_compressed.png'
    fig_fil = fig_fil.format(','.join(settings),stat)

    # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 25})
    else:
        plt.rcParams.update({'font.size': 100})

    if 'compressed' in fig_fil:
        fig,axs = plt.subplots(3,2,figsize=(31,33))
    else:
        fig,axs = plt.subplots(3,2,figsize=(140,133))
    #fig.suptitle('Hourly Errors and Values for each GRU', fontsize=40)
    fig.subplots_adjust(hspace=0.24) # Adjust the bottom margin, vertical space, and horizontal space

    # Specify variables of interest
    plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff']
    for i,(var,plt_t,leg_t,leg_t0,leg_tm) in enumerate(zip(plot_vars,plt_titl,leg_titl,leg_titl0,leg_titlm)): 
        run_loop(i,var,plt_t,leg_t,leg_t0,leg_tm)

    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)

if do_balance:
# Specify variables of interest
    plot_vars = ['balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceCasNrg','wallClockTime']
    comp_vars = ['balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','numberFluxCalc']
 
    plt_titl =  ['(a) Vegetation Balance','(b) Snow Balance','(c) Soil Balance', '(d) Canopy Air Space and Aquifer Balance', '(f) Wall Clock Time']
    leg_titl = ['$W~m^{-3}$'] * 4 + ['$s$']
    leg_titl0 =['$kg~m^{-2}~s^{-1}$'] * 4 + ['$num$']

    fig_fil = 'balance_scat_{}_compressed.png'
    if do_rel: fig_fil = 'balance_scat_{}_rel_compressed.png'
    fig_fil = fig_fil.format(stat)

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

    for i,(var,comp,leg_t,leg_t0,plt_t) in enumerate(zip(plot_vars,comp_vars,leg_titl,leg_titl0,plt_titl)): 
        run_loopb(i,var,comp,leg_t,leg_t0,plt_t)

    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
