# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python hist_per_GRU.py [stat]
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

do_rel = True # true is plot relative to the benchmark simulation
do_hist = False # true is plot histogram instead of CDF
run_local = True # true is run on local machine, false is run on cluster
fixed_Mass_units = False # true is convert mass balance units to kg m-2 s-1, if ran new code with depth in calculation

if run_local: 
    stat = 'rmnz'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')
    

#method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument
#plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6'] #maybe make this an argument
#method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument
#plt_name=['BE1','BE16','BE32','SUNDIALS'] #maybe make this an argument
method_name=['be1','be1cm','be1en','sundials_1en6cm'] 
plt_name=['BE1 common','BE1 temp','BE1 mixed','SUNDIALS temp']
method_name2=method_name+['sundials_1en8cm']
plt_name2=plt_name+['reference solution']

num_bins = 1000

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Define the power transformation function
def power_transform(x):
    return x ** 0.5  # Adjust the exponent as needed

# Simulation statistics file locations
use_vars = []
rep = [] # mark the repeats
use_vars = [1]
rep = [0] # mark the repeats
#use_vars = [0,1,2,3,4]
#rep = [0,0,0,0,0] # mark the repeats
settings0= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
settings = [settings0[i] for i in use_vars]

#use_vars2 = [4,4,5,5,6,6,7,7]
#use_vars2 = [0,0,1,1,2,2,3,3]
#rep2 = [1,2,1,2,1,2,1,2] # mark the repeats
use_vars2 = [0,0,1,1,2,2]
rep2 = [1,2,1,2,1,2] # mark the repeats
use_vars2 = [8,3,3]
rep2 = [0,1,2] # mark the repeats
#use_vars2 = [8]
#rep2 = [0] # mark the repeats
settings20= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']
settings2 = [settings20[i] for i in use_vars2]

viz_fil = method_name.copy()
viz_fl2 = method_name2.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings0))
for i, m in enumerate(method_name2):
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['balance']))

# Specify variables of interest
plot_vars = settings.copy()
plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','average routed runoff','wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$s$']
plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
leg_titl = [leg_titl[i] for i in use_vars]

plot_vars2 = settings2.copy()
plt_titl2 = ['canopy air space energy balance','vegetation energy balance','snow energy balance','soil energy balance','vegetation mass balance','snow mass balance','soil mass balance','aquifer mass balance', 'wall clock time']
leg_titl2 = ['$W~m^{-3}$'] * 4 + ['$kg~m^{-2}~s^{-1}$'] * 4 + ['$s$']
if fixed_Mass_units: leg_titl2 = ['$W~m^{-3}$'] * 4 + ['s^{-1}$'] * 3 + ['m~s^{-1}$'] + ['$s$']
plt_titl2 = [f"({chr(97+n + len(use_vars))}) {plt_titl2[i]}" for n,i in enumerate(use_vars2)]
leg_titl2 = [leg_titl2[i] for i in use_vars2]

if do_hist:
    fig_fil = 'Hrly_diff_hist_{}_{}_zoom_compressed.png'
    if do_rel: fig_fil = 'Hrly_diff_hist_{}_{}_zoom_rel_compressed.png'
else:
    fig_fil = 'Hrly_diff_cdf_{}_{}_zoom_compressed.png'
    if do_rel: fig_fil = 'Hrly_diff_cdf_{}_{}_zoom_rel_compressed.png'
fig_fil = fig_fil.format(','.join(settings),stat)

if stat == 'rmse' or stat=='rmnz':
    stat2 = 'mean'
    maxes = [2,15,250,0.08,200,10e-3]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,10e-3]
if stat == 'maxe':
    stat2 = 'amax'
    maxes = [15,25,0.8,2,0.3,2.0]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,2.0]
if stat == 'kgem':
    stat2 = 'mean'
    maxes = [0.9,0.9,0.9,0.9,0.9,10e-3]
maxes = [maxes[i] for i in use_vars]

if stat2 == 'mean':
    maxes2 = [1e-1,1e1,1e1,1e1]+[1e-7,1e-7,1e-7,1e-9] + [20e-3]
if stat2 == 'amax':
    maxes2 = [1e1,1e3,1e3,1e3]+[1e-5,1e-5,1e-5,1e-7] + [2.0]
maxes2 = [maxes2[i] for i in use_vars2]
for i in range(len(maxes2)):
    if rep2[i]==2: maxes2[i] = maxes2[i]*1e2 #clunky way to increase the range for the second repeat

summa = {}
summa1 = {}
if len(use_vars)>0:
    for i, m in enumerate(method_name):
        # Get the aggregated statistics of SUMMA simulations
        summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
if len(use_vars2)>0:
    for i, m in enumerate(method_name2):
        summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    
##Figure

plt.rcParams['xtick.color'] = 'black'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['ytick.major.width'] = 2

if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 27})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(4,2,figsize=(35,52))
else:
    fig,axs = plt.subplots(4,2,figsize=(140,160))
fig.subplots_adjust(hspace=0.33, wspace=0.17) # Adjust the bottom margin, vertical space, and horizontal space
#fig.suptitle('Histograms of Hourly Statistics for each GRU', fontsize=40,y=1.0)
    
def run_loop(i,var,mx,rep):
    r = i//2
    c = i-r*2
    stat0 = stat
    if rep == 1: stat0 = 'rmnz'
    if rep == 2: stat0 = 'maxe'
    if stat == 'rmse' or stat == 'kgem': 
        if var == 'wallClockTime': stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz':
        if var == 'wallClockTime': stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe': 
        if var == 'wallClockTime': stat0 = 'amax'
        statr = 'amax_ben'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx
    else:
        mx = 0.0
        mn = 1.0
        if do_rel: s_rel = summa[method_name[0]][var].sel(stat=statr)
        for m in method_name:
            s = summa[m][var].sel(stat=stat0)
            if do_rel and var != 'wallClockTime': s = s/s_rel
            if stat == 'maxe': s = np.fabs(s) # make absolute value norm
            mx = max(s.max(),mx)
            if stat == 'kgem': mn = min(s.min(),mn)

    # Data
    if do_rel: s_rel = summa[method_name[0]][var].sel(stat=statr)
    for m in method_name:
        s = summa[m][var].sel(stat=stat0)
        if do_rel and var != 'wallClockTime': s = s/s_rel

        if var == 'scalarTotalET' and not do_rel:
            if stat =='rmse' or stat =='rmnz' : s = s*31557600 # make annual total
            if stat =='maxe': s = s*3600 # make hourly max
        if var == 'averageRoutedRunoff' and not do_rel:
            if stat =='rmse' or stat =='rmnz' : s = s*31557600*1000 # make annual total
            if stat =='maxe': s = s*3600*1000 # make hourly max           
        if stat == 'maxe': s = np.fabs(s) # make absolute value norm
        range = (0,mx)
        if stat=='kgem' and var!='wallClockTime' : range = (mn,1)
        if do_hist: 
            np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=range)
        else: #cdf
            sorted_data = np.sort(np.fabs(s))
            valid_data = sorted_data[~np.isnan(sorted_data)]
            yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
            axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=2.0)
            axs[r,c].set_xlim(range)  # Replace xmin and xmax with the desired limits


    if stat0 == 'rmse': stat_word = 'RMSE'
    if stat0 == 'rmnz': stat_word = 'RMSE' # no 0s'
    if stat0 == 'maxe': stat_word = 'max abs error'
    if stat0 == 'kgem': stat_word = 'KGE"'
    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'mnnz': stat_word = 'mean' # no 0s'
    if stat0 == 'amax': stat_word = 'max'
    
    if statr == 'mean_ben': statr_word = 'mean'
    if statr == 'mnnz_ben': statr_word = 'mean' # no 0s'
    if statr == 'amax_ben': statr_word = 'max'
    
    axs[r,c].legend(plt_name)
    axs[r,c].set_title(plt_titl[i])
    if rep>0: axs[r,c].set_title(plt_titl[i] + ' '+ stat_word)
    if stat == 'rmse' or stat == 'rmnz' or stat == 'maxe': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl[i]))
    if stat == 'kgem': axs[r,c].set_xlabel(stat_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    if do_hist: 
        axs[r,c].set_ylabel('GRU count')
        if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
 
    else:
        axs[r,c].set_ylabel('cumulative distribution')
        if(c==1): axs[r, c].set_ylabel('')
        axs[r,c].set_ylim([0.0, 1.0])
        axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
        if var=='scalarTotalSoilWat' or var=='wallClockTime': # Rotate x-axis labels for axs[2, 1] subplot
            axs[r, c].tick_params(axis='x', rotation=45)

def run_loopb(i,var,mx,rep):
    r = (i+len(use_vars))//2
    c = (i+len(use_vars))-r*2
    stat0 = stat2
    if rep == 1: stat0 = 'mean'
    if rep == 2: stat0 = 'amax'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx*1e-9
        if any(substring in var for substring in ['VegNrg', 'SnowNrg', 'SoilNrg']):
            mn = mx*1e-9
        if var=='wallClockTime': mn = 0.0
        if fixed_Mass_units and 'Mass' in var: # /density for mass balance
            mn = mn/1000
            mx = mx/1000
    else:
        mx = 0.0
        mn = 1.0
        for m in method_name2:
            # Get the statistics, remove 9999 (should be nan, but just in case)
            s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
            mx = max(s.max(),mx)
            mn = min(s.min(),mn)

    # Data
    for m in method_name2:
        s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
        if fixed_Mass_units and 'Mass' in var: s = s/1000 # /density for mass balance

        range = (mn,mx)
        if do_hist: 
            np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=range)
        else: #cdf
            sorted_data = np.sort(np.fabs(s))
            valid_data = sorted_data[~np.isnan(sorted_data)]
            yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
            axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=2.0)
            axs[r,c].set_xlim(range)  # Replace xmin and xmax with the desired limits


    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'amax': stat_word = 'max'

    axs[r,c].legend(plt_name2)
    axs[r,c].set_title(plt_titl2[i])
    if rep>0: axs[r,c].set_title(plt_titl2[i] + ' '+ stat_word)
    axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl2[i]))   

    if do_hist: 
        axs[r,c].set_ylabel('GRU count')
        if(c==1): axs[r, c].set_ylabel('')
        if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
 
    else:
        axs[r,c].set_ylabel('cumulative distribution')
        if(c==1): axs[r, c].set_ylabel('')
        axs[r,c].set_ylim([0.0, 1.0])
        axs[r,c].set_xscale('log') #log x axis
        if var=='wallClockTime': 
            axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
            axs[r, c].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot

if len(use_vars) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars,maxes,rep)): 
        run_loop(i,var,mx,rep)
if len(use_vars2) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars2,maxes2,rep2)): 
        run_loopb(i,var,mx,rep)

# Remove the extra subplots
if (len(plot_vars)+len(plot_vars2)) < 8:
    for i in range((len(plot_vars)+len(plot_vars2)),8):
        r = i//2
        c = i-r*2
        fig.delaxes(axs[r, c])

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
