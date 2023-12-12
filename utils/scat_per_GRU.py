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
use_eff = False # use efficiency in wall clock time, still need files for the node number
do_rel = True # plot relative to the benchmark simulation

testing = False
if testing: 
    stat = 'maxe'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1','sundials_1en6'] #maybe make this an argument
    plt_name=['BE1','IDAe-6'] #maybe make this an argument
else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    #method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument
    #plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6'] #maybe make this an argument
    method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument
    plt_name=['BE1','BE16','BE32','IDAe-6'] #maybe make this an argument

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
eff_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
    eff_fil[i] = 'eff_' + m + '.txt'

# Specify variables of interest
plot_vars = settings.copy()
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$num$']
leg_titl0 = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$s$']
leg_titlm= ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~h^{-1}$','$kg~m^{-2}$','$mm~h^{-1}$','$num$']

fig_fil = 'Hrly_diff_scat_{}_{}_compressed.png'
if do_rel: fig_fil = 'Hrly_diff_scat_{}_{}_rel_compressed.png'
fig_fil = fig_fil.format(','.join(settings),stat)

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

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(3,2,figsize=(31,33))
else:
    fig,axs = plt.subplots(3,2,figsize=(140,133))
#fig.suptitle('Hourly Errors and Values for each GRU', fontsize=40)
fig.subplots_adjust(hspace=0.24) # Adjust the bottom margin, vertical space, and horizontal space

    
def run_loop(i,var):
    r = i//2
    c = i-r*2
    if stat == 'rmse' or stat == 'kgem': 
        stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz':
        stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe': 
        stat0 = 'amax'
        statr = 'amax_ben'

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

        if var == 'wallClockTime':
            batch = np.floor(np.arange(len(s.indexes['hru'])) /nbatch_hrus)
            #basin_num = np.arange(len(s.indexes['hru'])) % nbatch_hrus #not currently using
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
            node_batch = np.array([node[b] for b in batch]) #not currently using
            # Multiply the s values by efficiency
            if use_eff: s = s*eff_batch
            axs[r,c].scatter(x=node_batch,y=s.sel(stat=stat0),s=1,zorder=0,label=m)
            stat_word = 'Node number'
        else:
            axs[r,c].scatter(x=np.fabs(s.sel(stat=stat).values),y=s.sel(stat=stat0).values,s=1,zorder=0,label=m)        
            if stat == 'rmse': stat_word = 'RMSE'
            if stat == 'rmnz': stat_word = 'RMSE' # no 0s'
            if stat == 'maxe': stat_word = 'max abs error'
            if stat == 'kgem': stat_word = 'KGE"'
 
    if stat0 == 'mean': stat0_word = 'mean'
    if stat0 == 'mnnz': stat0_word = 'mean' # no 0s'
    if stat0 == 'amax': stat0_word = 'max'
 
    lgnd = axs[r,c].legend(plt_name)
    for j, m in enumerate(method_name):
       lgnd.legendHandles[j]._sizes = [80]
    axs[r,c].set_title(plt_titl[i])
    if stat == 'rmse' or stat == 'rmnz': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl[i]))
    if stat == 'maxe': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titlm[i]))   
    if stat == 'kgem': axs[r,c].set_xlabel(stat_word)
    #if do_rel and var!='wallClockTime': axs[r,c].set_xlabel(stat_word + ' rel to bench ' + stat0_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    axs[r,c].set_ylabel(stat0_word + ' [{}]'.format(leg_titl0[i]))
    #if do_rel and var!='wallClockTime': axs[r,c].set_ylabel(stat0_word + ' rel to bench ' + stat0_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_ylabel('relative '+ stat0_word)

for i,var in enumerate(plot_vars): 
    run_loop(i,var)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
