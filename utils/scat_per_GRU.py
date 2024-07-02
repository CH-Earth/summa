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
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap

viz_dir = Path('/home/avanb/scratch/statistics')
nbatch_hrus = 518 # number of HRUs per batch
do_rel = False # plot relative to the benchmark simulation
do_heat = True # plot heatmaps instead of scatterplots
rainbow_cmap = plt.cm.get_cmap('rainbow', 256)  # Get the rainbow colormap
rainbow_colors = rainbow_cmap(np.linspace(0, 1, 256))
rainbow_colors_with_white = np.vstack((rainbow_colors, [1, 1, 1, 1]))
custom_cmap = ListedColormap(rainbow_colors_with_white, name='rainbow_white')
custom_cmap.set_under('white')  # Ensure that values under the lower bound are white

# which statistics to plot, can do both
do_vars = True
do_balance = True

testing = False
if testing: 
    stat = 'rmnz'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
    method_name=['be1en']
    plt_name=['BE1 mixed']
    method_name2=method_name
    plt_name2=plt_name

else:
    import sys
    # The first input argument specifies the run where the files are
    stat = sys.argv[1]
    #method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument
    #plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6'] #maybe make this an argument
    #method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument
    #plt_name=['BE1','BE16','BE32','IDAe-6'] #maybe make this an argument
    method_name=['be1','be1cm','be1en','sundials_1en6cm'] 
    plt_name=['BE1 common','BE1 temp','BE1 mixed','SUNDIALS temp']
    method_name2=method_name+['sundials_1en8cm']
    plt_name2=plt_name+['reference solution']

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_fil = method_name.copy()
viz_fl2 = method_name2.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
for i, m in enumerate(method_name2):
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['balance']))

summa = {}
summa1 = {}
if do_vars:
    for i, m in enumerate(method_name):
        # Get the aggregated statistics of SUMMA simulations
        summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
if do_balance:
    for i, m in enumerate(method_name2):
        summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    
def run_loop(i,var,lx,ly,plt_t,leg_t,leg_t0,leg_tm):

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

    if stat == 'mnnz':
        stat0 = 'mnnz_ben'
        stat0_word = 'mean - mean reference' # no 0s'
        statr = 'mnnz_ben'
        stat_word = 'mean' # no 0s'

    if stat == 'mean':
        stat0 = 'mean_ben'
        stat0_word = 'mean - mean reference' # no 0s'
        statr = 'mean_ben'
        stat_word = 'mean' # no 0s'

    if stat == 'amax':
        stat0 = 'amax_ben'
        stat0_word = 'max abs - max abs reference'
        statr = 'amax_ben'
        stat_word = 'max abs error'

    # Data
    if do_rel: s_rel = summa[method_name[0]][var].sel(stat=statr)
    for j, m in enumerate(method_name):
        s = summa[m][var].sel(stat=[stat,stat0])
        if var == 'scalarTotalET':
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz' or stat=='mean': 
                s = s*31557600 # make annual total
                if do_rel: s_rel = s_rel*31557600
            if stat =='maxe' or stat=='amax': 
                s = s*3600 # make hourly max
                if do_rel: s_rel = s_rel*3600
        if var == 'averageRoutedRunoff':
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz'or stat=='mean': 
                s = s*31557600*1000 # make annual total
                if do_rel: s_rel = s_rel*31557600*1000
            if stat =='maxe' or stat=='amax': 
                s = s*3600*1000 # make hourly max
                if do_rel: s_rel = s_rel*3600*1000
        if stat == 'maxe': s.loc[dict(stat='maxe')] = np.fabs(s.loc[dict(stat='maxe')]) # make absolute value norm

        if do_rel and var != 'wallClockTime': 
            s.loc[dict(stat=stat)] = s.loc[dict(stat=stat)]/s_rel
            stat_word = 'relative '+ stat_word

        if do_heat and j==0:
            x_points = []
            y_points = []
            if stat=='mnnz' or stat=='mean' or stat=='amax': 
                x = s.sel(stat=stat).values
                y = s.sel(stat=stat).values-s.sel(stat=stat0).values
            else:
                x = np.fabs(s.sel(stat=stat).values)
                y = s.sel(stat=stat0).values
            if lx:
                x = np.where(x > 0, np.log10(x), np.nan)
                stat_word = 'log10 ' + stat_word
            if ly:
                if stat!='mnnz' and stat!='mean': 
                    if var == 'scalarTotalET':
                        y = np.where(-y > 0, np.log10(-y), np.nan)
                        stat0_word = 'log10 negative ' + stat0_word
                    else:
                        y = np.where(y > 0, np.log10(y), np.nan)
                        stat0_word = 'log10 ' + stat0_word
            x_points.extend(x)
            y_points.extend(y)
            x_points = np.array(x_points)  # Convert lists to numpy arrays
            y_points = np.array(y_points)

            # Ensure no NaNs or infs before proceeding
            is_valid = ~np.isnan(x_points) & ~np.isnan(y_points) & ~np.isinf(x_points) & ~np.isinf(y_points)
            x_points = x_points[is_valid]
            y_points = y_points[is_valid]
            if not is_valid.any():
                print('no valid values')
                continue
            numbin = int(np.sqrt(x_points.size)/2) + 1
            print(var,'numbin', numbin)

            # Define the bins for the histogram and calculate
            x_edges = np.linspace(x_points.min(), x_points.max(), num=numbin)
            y_edges = np.linspace(y_points.min(), y_points.max(), num=numbin)
            zi, _, _ = np.histogram2d(x_points, y_points, bins=[x_edges, y_edges])

            # Calculate bin centers from edges for X and Y
            x_centers = (x_edges[:-1] + x_edges[1:]) / 2
            y_centers = (y_edges[:-1] + y_edges[1:]) / 2
            X, Y = np.meshgrid(x_centers, y_centers)

            # Adjust the pcolormesh call to use the centers and compatible shading
            norm = LogNorm(vmin=np.min(zi[zi > 0]), vmax=np.max(zi))
            mesh = axs[r, c].pcolormesh(X, Y, zi.T, shading='gouraud', cmap=custom_cmap, zorder=0,norm=norm)
            fig.colorbar(mesh, ax=axs[r, c], label='GRU count')

        elif not do_heat:
            if stat=='mnnz' or stat=='mean' or stat=='amax': 
                axs[r,c].scatter(x=s.sel(stat=stat).values,y=s.sel(stat=stat).values-s.sel(stat=stat0).values,s=1,zorder=0,label=m)       
            else: 
                axs[r,c].scatter(x=np.fabs(s.sel(stat=stat).values),y=s.sel(stat=stat0).values,s=1,zorder=0,label=m)        

    if do_heat:
        axs[r,c].set_title(plt_t + ' '+ plt_name[0] + ' heatmap')
    elif not(do_heat):
        lgnd = axs[r,c].legend(plt_name)
        for j, m in enumerate(plt_name):
           lgnd.legendHandles[j]._sizes = [80]
        axs[r,c].set_title(plt_t)
    if do_rel: 
        axs[r,c].set_xlabel(stat_word)  
    else:
        if stat == 'rmse' or stat == 'rmnz' or stat=='mnnz' or stat=='mean': 
            axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_t))
        if stat == 'maxe' or stat=='amax': 
            axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_tm)) 
    if stat == 'kgem': axs[r,c].set_xlabel(stat_word)
    axs[r,c].set_ylabel(stat0_word + ' [{}]'.format(leg_t0))

def run_loopb(i,var,comp,lx,ly,leg_t,leg_t0,plt_t):
    r = i//2
    c = i-r*2

    if stat == 'rmse' or stat == 'kgem': 
        stat0 = 'mean'
        wordx = ' mean'
    if stat == 'rmnz':
        stat0 = 'mean'
        wordx = ' mean' # no 0s'
    if stat == 'maxe': 
        stat0 = 'amax'
        wordx = ' max'
    wordy = np.copy(wordx)

    # Data
    for j, m in enumerate(method_name2):
        # Get the statistics, remove 9999 (should be nan, but just in case)
        s0 = np.fabs(summa1[m][comp].sel(stat=stat0)).where(lambda x: x != 9999)
        s = np.fabs(summa1[m][var].sel(stat=stat0)).where(lambda x: x != 9999)

        if do_heat and j==0:
            x_points = []
            y_points = []
            x = s.values
            y = s0.values
            x_points.extend(x)
            y_points.extend(y)
            x_points = np.array(x_points)  # Convert lists to numpy arrays
            y_points = np.array(y_points)
            if lx:
                x = np.log10(x[x > 0])
                wordx = 'log10 ' + wordx
            if ly:
                y = np.log10(y[y > 0])
                wordy = 'log10 ' + wordy

            x_points.extend(x)
            y_points.extend(y)
            x_points = np.array(x_points)  # Convert lists to numpy arrays
            y_points = np.array(y_points)

            # Ensure no NaNs or infs before proceeding
            is_valid = ~np.isnan(x_points) & ~np.isnan(y_points) & ~np.isinf(x_points) & ~np.isinf(y_points)
            x_points = x_points[is_valid]
            y_points = y_points[is_valid]
            numbin = int(np.sqrt(x_points.size)/2) + 1
            print(var,'numbin', numbin)

            # Define the bins for the histogram and calculate
            x_edges = np.linspace(x_points.min(), x_points.max(), num=numbin)
            y_edges = np.linspace(y_points.min(), y_points.max(), num=numbin)
            zi, _, _ = np.histogram2d(x_points, y_points, bins=[x_edges, y_edges])

            # Calculate bin centers from edges for X and Y
            x_centers = (x_edges[:-1] + x_edges[1:]) / 2
            y_centers = (y_edges[:-1] + y_edges[1:]) / 2
            X, Y = np.meshgrid(x_centers, y_centers)

            # Adjust the pcolormesh call to use the centers and compatible shading
            norm = LogNorm(vmin=np.min(zi[zi > 0]), vmax=np.max(zi))
            mesh = axs[r, c].pcolormesh(X, Y, zi.T, shading='gouraud', cmap=custom_cmap, zorder=0,norm=norm)
            fig.colorbar(mesh, ax=axs[r, c], label='GRU count')

        elif not do_heat:
            axs[r,c].scatter(x=s.values,y=s0.values,s=10,zorder=0,label=m)        

    if comp == 'numberFluxCalc':
        stat0_word = 'number flux calculations'
        stat_word = 'wall clock time'
    else:
        stat0_word = 'balance abs value'
        stat_word = 'balance abs value'

    if do_heat:
        axs[r,c].set_title(plt_t + ' '+ plt_name[0] + ' heatmap')
    elif not(do_heat):
        lgnd = axs[r,c].legend(plt_name2)
        for j, m in enumerate(plt_name2):
           lgnd.legendHandles[j]._sizes = [80]
        axs[r,c].set_title(plt_t)
        axs[r,c].set_xscale('log')
        if comp != 'numberFluxCalc': axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(stat_word  + wordx + ' [{}]'.format(leg_t))
    axs[r,c].set_ylabel(stat0_word + wordy + ' [{}]'.format(leg_t0))

plt.rcParams['xtick.color'] = 'black'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['ytick.major.width'] = 2
if do_vars:
    # Specify variables of interest
    use_vars = [1,2,4]
    #logx = [0,0,0] # log scale x axis
    #logy = [0,0,1] # log scale y axis
    use_vars = [0,1,2,3,4]
    logx = np.ones(len(use_vars)) # log scale x axis
    logy = np.ones(len(use_vars)) # log scale y axis

    plot_vars = settings
    plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','average routed runoff']
    leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$mm~y^{-1}$']
    leg_titl0 = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$mm~y^{-1}$']
    leg_titlm= ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~h^{-1}$','$kg~m^{-2}$','$mm~h^{-1}$','$mm~h^{-1}$']

    plot_vars = [plot_vars[i] for i in use_vars]
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
    leg_titl = [leg_titl[i] for i in use_vars]
    leg_titl0 = [leg_titl0[i] for i in use_vars]
    leg_titlm = [leg_titlm[i] for i in use_vars]

    fig_fil = 'Hrly_diff_scat_{}_{}'
    if do_rel: fig_fil = fig_fil + '_rel'
    if do_heat: fig_fil = method_name[0]+fig_fil + '_heat'
    fig_fil = fig_fil + '_compressed.png'
    fig_fil = fig_fil.format(','.join(settings),stat)

    # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 27})
    else:
        plt.rcParams.update({'font.size': 100})

    if 'compressed' in fig_fil:
        fig,axs = plt.subplots(3,2,figsize=(35,38))
    else:
        fig,axs = plt.subplots(3,2,figsize=(140,133))
    #fig.suptitle('Hourly Errors and Values for each GRU', fontsize=40)
    fig.subplots_adjust(hspace=0.33, wspace=0.17) # Adjust the bottom margin, vertical space, and horizontal space

    for i,(var,lx,ly,plt_t,leg_t,leg_t0,leg_tm) in enumerate(zip(plot_vars,logx,logy,plt_titl,leg_titl,leg_titl0,leg_titlm)): 
        run_loop(i,var,lx,ly,plt_t,leg_t,leg_t0,leg_tm)

    # Remove the extra subplots
    if (len(plot_vars)) < 6:
        for i in range(len(plot_vars),6):
            r = i//2
            c = i-r*2
            fig.delaxes(axs[r, c])

    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)

if do_balance:
# Specify variables of interest
    use_vars = [0,1,2,3]
    logx = np.ones(len(use_vars)) # log scale x axis
    logy = np.ones(len(use_vars)) # log scale y axis

    plot_vars = ['balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceCasNrg','wallClockTime']
    comp_vars = ['balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','numberFluxCalc']
    plt_titl = ['vegetation balance','snow balance','soil balance', 'canopy air space and aquifer balance', 'wall clock time']
    leg_titl = ['$W~m^{-3}$'] * 4 + ['$s$']
    leg_titl0 =['$kg~m^{-2}~s^{-1}$'] * 4 + ['$num$']

    plot_vars = [plot_vars[i] for i in use_vars]
    comp_vars = [comp_vars[i] for i in use_vars]
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
    leg_titl = [leg_titl[i] for i in use_vars]
    leg_titl0 = [leg_titl0[i] for i in use_vars]

    fig_fil = 'balance_scat_{}'
    if do_heat: fig_fil = method_name[0]+fig_fil + '_heat'
    fig_fil = fig_fil + '_compressed.png'
    fig_fil = fig_fil.format(stat)

    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 27})
    else:
        plt.rcParams.update({'font.size': 100})

    if 'compressed' in fig_fil:
        fig,axs = plt.subplots(3,2,figsize=(35,38))
    else:
        fig,axs = plt.subplots(3,2,figsize=(140,160))
    fig.subplots_adjust(hspace=0.33, wspace=0.17) # Adjust the bottom margin, vertical space, and horizontal space
    #fig.suptitle('Scatterplot of Hourly Statistics for each GRU', fontsize=40,y=1.0)

    for i,(var,comp,lx,ly,leg_t,leg_t0,plt_t) in enumerate(zip(plot_vars,comp_vars,logx,logy,leg_titl,leg_titl0,plt_titl)): 
        run_loopb(i,var,comp,lx,ly,leg_t,leg_t0,plt_t)

    if (len(plot_vars)) < 6:
        for i in range(len(plot_vars),6):
            r = i//2
            c = i-r*2
            fig.delaxes(axs[r, c])
    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
