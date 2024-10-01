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

do_rel = False # true is plot relative to the benchmark simulation
do_heat = False # true is plot heatmaps instead of scatterplots
run_local = True # true is run on local machine, false is run on cluster
inferno_col= True # Set to True if want to match geographic plots, False if want rainbow colormap (does not matter if do_heat is False)
fixed_Mass_units = False # true is convert mass balance units to kg m-2 s-1, if ran new code with depth in calculation

# which statistics to plot, can do both
do_vars = False
do_balance = True

if run_local: 
    stat = 'mnnz'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')


#method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6']
#plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6']
method_name=['be1','be16','be32','sundials_1en6']
plt_name=['BE1','BE16','BE32','SUNDIALS']
method_name=['be1','be1cm','be1en','sundials_1en6cm','sundials_1en6en'] 
plt_name=['BE1 common','BE1 temp','BE1 mixed','SUNDIALS temp','SUNDIALS enth']
method_name2=method_name+['sundials_1en8cm']
plt_name2=plt_name+['reference solution']

if inferno_col:
    custom_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    custom_cmap.set_bad(color='white') #nan color white
else: # use rainbow colormap, I think looks better
    rainbow_cmap = plt.cm.get_cmap('rainbow', 256)  # Get the rainbow colormap
    rainbow_colors = rainbow_cmap(np.linspace(0, 1, 256))
    rainbow_colors_with_white = np.vstack((rainbow_colors, [1, 1, 1, 1]))
    custom_cmap = ListedColormap(rainbow_colors_with_white, name='rainbow_white')
    custom_cmap.set_under('white')  # Ensure that values under the lower bound are white

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
    hru_size = summa[m].sizes['hru']
if do_balance:
    for i, m in enumerate(method_name2):
        summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    hru_size = summa1[m].sizes['hru']

numbin = int(np.sqrt(hru_size/10))
maxcolor = numbin**2/75
do_clip = True # choose if want the heat values clipped (True) or plotted as white if over maxcolor (False)
    
def run_loop(i,var,lx,ly,plt_t,leg_t,leg_t0,leg_tm,rep,mx):
    r = i//ncol
    c = i-r*ncol
    global method_name  # Declare method_name as global to modify its global value
    global plt_name  # Declare plt_name as global to modify its global value
    method_name = np.copy(method_name0)
    plt_name = np.copy(plt_name0)

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
    # make the axes the same
    do_same = False
    if do_heat and len(method_name)>1:
            do_same = True
            mxx = 0.0
            mnx = 1.0
            for m in method_name:
                # Get the statistics, remove 9999 (should be nan, but just in case)
                s0 = summa[m][var].sel(stat=stat).where(lambda x: x != 9999)
                if do_rel: s0=s0/s_rel
                if stat=='mnnz' or stat=='mean' or stat=='amax':
                    s = s0
                else:
                    s = np.fabs(s0)
                mxx = max(s.max(),mxx)
                mnx = min(s.min(),mnx)
            mxy = 0.0
            mny = 1.0
            for m in method_name:
                # Get the statistics, remove 9999 (should be nan, but just in case)
                s0 = summa[m][var].sel(stat=[stat,stat0]).where(lambda x: x != 9999)
                if stat=='mnnz' or stat=='mean' or stat=='amax':
                    s = s0.sel(stat=stat) - s0.sel(stat=stat0)
                else:
                    s = s0.sel(stat=stat0)
                mxy = max(s.max(),mxy)
                mny = min(s.min(),mny)
            method_name = [method_name0[rep]]
            plt_name = [plt_name0[rep]]
    else: # only one method
        mxx = 0.0
        mnx = 1.0
        mxy = 0.0
        mny = 1.0
    for j, m in enumerate(method_name):
        s = summa[m][var].sel(stat=[stat,stat0])
        if var == 'scalarTotalET':
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz' or stat=='mean': 
                s = s*31557600 # make annual total
                if do_rel: 
                    s_rel = s_rel*31557600
                else:
                    mnx = mnx*31557600
                    mxx = mxx*31557600
                mny = mny*31557600
                mxy = mxy*31557600
            if stat =='maxe' or stat=='amax': 
                s = s*3600 # make hourly max
                if do_rel: 
                    s_rel = s_rel*3600
                else:
                    mnx = mnx*3600
                    mxx = mxx*3600
                mny = mny*3600
                mxy = mxy*3600
        if var == 'averageRoutedRunoff':
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz'or stat=='mean': 
                s = s*31557600*1000 # make annual total
                if do_rel: 
                    s_rel = s_rel*31557600*1000
                else:
                    mnx = mnx*31557600*1000
                    mxx = mxx*31557600*1000
                mny = mny*31557600*1000
                mxy = mxy*31557600*1000
            if stat =='maxe' or stat=='amax': 
                s = s*3600*1000 # make hourly max
                if do_rel: 
                    s_rel = s_rel*3600*1000
                else:
                    mnx = mnx*3600*1000
                    mxx = mxx*3600*1000
                mny = mny*3600*1000
                mxy = mxy*3600*1000
        if stat == 'maxe': s.loc[dict(stat='maxe')] = np.fabs(s.loc[dict(stat='maxe')]) # make absolute value norm

        if do_rel and var != 'wallClockTime': 
            s.loc[dict(stat=stat)] = s.loc[dict(stat=stat)]/s_rel
            stat_word = 'relative '+ stat_word

        if do_heat:
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
                mnx = np.log10(np.where(mnx <= 0, 1e-30, mnx))
                mxx = np.log10(np.where(mxx <= 0, 1e-30, mxx))
            if ly:
                if stat!='mnnz' and stat!='mean': 
                    if var == 'scalarTotalET':
                        y = np.where(-y > 0, np.log10(-y), np.nan)
                        mny = np.log10(np.where(-mny <= 0, 1e-30, -mny))
                        mxy = np.log10(np.where(-mxy <= 0, 1e-30, -mxy))
                    else:
                        y = np.where(y > 0, np.log10(y), np.nan)
                        mny = np.log10(np.where(mny <= 0, 1e-30, mny))
                        mxy = np.log10(np.where(mxy <= 0, 1e-30, mxy))
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

            # Define the bins for the histogram and calculate
            if not do_same:
                mnx = x_points.min()
                mxx = x_points.max()
                mny = y_points.min()
                mxy = y_points.max()
            else:
                if mx:
                    mnx = 0.0
                    mxx = mx
            x_edges = np.linspace(mnx,mxx, num=numbin)
            y_edges = np.linspace(mny,mxy, num=numbin)            
            zi, _, _ = np.histogram2d(x_points, y_points, bins=[x_edges, y_edges])
            zi = zi.T  # Transpose the histogram to match the pcolormesh orientation
            if do_clip: 
                zi_clipped = np.where(zi > 0.95*maxcolor, 0.95*maxcolor, zi) # Clip zi.T so that max values are maxcolor with a little buffer
            else:
                zi_clipped = zi # don't clip

            # Calculate bin centers from edges for X and Y
            x_centers = (x_edges[:-1] + x_edges[1:]) / 2
            y_centers = (y_edges[:-1] + y_edges[1:]) / 2
            X, Y = np.meshgrid(x_centers, y_centers)
            if lx: X = 10**X
            if ly: 
                if stat!='mnnz' and stat!='mean': 
                    if var == 'scalarTotalET':
                        Y = -10**Y
                    else:
                        Y = 10**Y

            # Adjust the pcolormesh call to use the centers and compatible shading
            norm = LogNorm(vmin=1, vmax=maxcolor)
            mesh = axs[r, c].pcolormesh(X, Y, zi_clipped, shading='gouraud', cmap=custom_cmap, zorder=0,norm=norm)
            if r==1 and c==len(method_name)-1: fig.colorbar(mesh, ax=axs.ravel().tolist(), label='GRU count',aspect=20/3*nrow)

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
    if lx: axs[r,c].set_xscale('log')
    if ly: axs[r,c].set_yscale('log')
    if do_rel: 
        axs[r,c].set_xlabel(stat_word)  
    else:
        if stat == 'rmse' or stat == 'rmnz' or stat=='mnnz' or stat=='mean': 
            axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_t))
        if stat == 'maxe' or stat=='amax': 
            axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_tm)) 
    if stat == 'kgem': axs[r,c].set_xlabel(stat_word)
    axs[r,c].set_ylabel(stat0_word + ' [{}]'.format(leg_t0))
    if do_heat and c>0: axs[r, c].set_ylabel('')

def run_loopb(i,var,comp,lx,ly,leg_t,leg_t0,plt_t,repy):
    r = i//ncol
    c = i-r*ncol
    global method_name2  # Declare method_name as global to modify its global value
    global plt_name2 # Declare plt_name as global to modify its global value
    method_name2 = np.copy(method_name20)
    plt_name2 = np.copy(plt_name20)

    if stat == 'rmse' or stat == 'kgem' or stat == 'mean': 
        stat0 = 'mean'
        wordx = ' mean'
    if stat == 'rmnz' or stat == 'mnnz':
        stat0 = 'mean'
        wordx = ' mean' # no 0s'
    if stat == 'maxe': 
        stat0 = 'amax'
        wordx = ' max'
    wordy = wordx

    # Data
    do_same = False
    if do_heat and len(method_name2)>1:
        if len(method_name2)>1:
            do_same = True
            mxx = 0.0
            mnx = 1.0
            for m in method_name2:
                # Get the statistics, remove 9999 (should be nan, but just in case)
                s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
                mxx = max(s.max(),mxx)
                mnx = min(s.min(),mnx)
            mxy = 0.0
            mny = 1.0
            for m in method_name2:
                # Get the statistics, remove 9999 (should be nan, but just in case)
                s = summa1[m][comp].sel(stat=stat0).where(lambda x: x != 9999)
                mxy = max(s.max(),mxy)
                mny = min(s.min(),mny)
            method_name2 = [method_name20[rep]]
            plt_name2 = [plt_name20[rep]]
    else: # only one method
        mxx = 0.0
        mnx = 1.0
        mxy = 0.0
        mny = 1.0
    for j, m in enumerate(method_name2):
        # Get the statistics, remove 9999 (should be nan, but just in case)
        s = np.fabs(summa1[m][var].sel(stat=stat0)).where(lambda x: x != 9999, np.nan)
        s0 = np.fabs(summa1[m][comp].sel(stat=stat0)).where(lambda x: x != 9999, np.nan)
        if fixed_Mass_units and 'Mass' in comp: s0 = s0/1000 # / density for mass balance

        if do_heat:
            x_points = []
            y_points = []
            x = s.values
            y = s0.values
            if lx:
                x = np.where(x > 0, np.log10(x), np.nan)
                mxx = np.log10(np.where(mxx <= 0, 1e-30, mxx))
                mnx = np.log10(np.where(mnx <= 0, 1e-30, mnx))
            if ly:
                y = np.where(y > 0, np.log10(y), np.nan)
                mxy = np.log10(np.where(mxy <= 0, 1e-30, mxy))
                mny = np.log10(np.where(mny <= 0, 1e-30, mny))
            x_points.extend(x)
            y_points.extend(y)
            x_points = np.array(x_points)  # Convert lists to numpy arrays
            y_points = np.array(y_points)

            # Ensure no NaNs or infs before proceeding
            is_valid = ~np.isnan(x_points) & ~np.isnan(y_points) & ~np.isinf(x_points) & ~np.isinf(y_points)
            x_points = x_points[is_valid]
            y_points = y_points[is_valid]

            # Define the bins for the histogram and calculate
            if not do_same:
                mnx = x_points.min()
                mxx = x_points.max()
                mny = y_points.min()
                mxy = y_points.max()
            x_edges = np.linspace(mnx,mxx, num=numbin)
            y_edges = np.linspace(mny,mxy, num=numbin)
            zi, _, _ = np.histogram2d(x_points, y_points, bins=[x_edges, y_edges])
            zi = zi.T  # Transpose the histogram to match the pcolormesh orientation
            if do_clip: 
                zi_clipped = np.where(zi > 0.95*maxcolor, 0.95*maxcolor, zi) # Clip zi.T so that max values are maxcolor with a little buffer
            else:
                zi_clipped = zi # don't clip
            # Calculate bin centers from edges for X and Y
            x_centers = (x_edges[:-1] + x_edges[1:]) / 2
            y_centers = (y_edges[:-1] + y_edges[1:]) / 2
            X, Y = np.meshgrid(x_centers, y_centers)
            if lx: X = 10**X
            if ly: Y = 10**Y

            # Adjust the pcolormesh call to use the centers and compatible shading
            norm = LogNorm(vmin=1, vmax=maxcolor)
            mesh = axs[r, c].pcolormesh(X, Y, zi_clipped, shading='gouraud', cmap=custom_cmap, zorder=0,norm=norm)
            if r==1 and c==len(method_name2)-1: fig.colorbar(mesh, ax=axs.ravel().tolist(), label='GRU count',aspect=20/3*nrow)
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
    if lx: axs[r,c].set_xscale('log')
    if ly: axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(stat_word  + wordx + ' [{}]'.format(leg_t))
    axs[r,c].set_ylabel(stat0_word + wordy + ' [{}]'.format(leg_t0))
    if do_heat and c>0: axs[r, c].set_ylabel('')

plt.rcParams['xtick.color'] = 'black'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['ytick.major.width'] = 2
if do_vars:
    # Specify variables of interest
    if do_heat:
        use_vars = [0,1,2,3,4]
        use_meth = [0,2]
        logx = np.zeros(len(use_vars)) # no log scale x axis
        logy = np.zeros(len(use_vars)) # log scale y axis
    else:
        use_vars = [0,1,2,3,4]
        use_meth = [0,1,2,3]
        logx = np.ones(len(use_vars)) # log scale x axis
        logy = np.ones(len(use_vars)) # log scale y axis

    rep = np.zeros(len(use_vars)) 
    if do_heat:
        use_vars = [val for val in use_vars for _ in range(len(use_meth))]
        logy = [val for val in logy for _ in range(len(use_meth))]
        logx = [val for val in logx for _ in range(len(use_meth))]
        rep = [int(val+_) for val in rep for _ in range(len(use_meth))]
             
    plot_vars = settings
    plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','average routed runoff']
    leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$mm~y^{-1}$']
    leg_titl0 = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$mm~y^{-1}$']
    leg_titlm= ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~h^{-1}$','$kg~m^{-2}$','$mm~h^{-1}$','$mm~h^{-1}$']

    #to zoom the heat x axis set these
    maxes = [0.0,0.0,0.0,0.0,0.0,0.0] #initialize
    if stat == 'rmse' or stat=='rmnz':
        maxes = [60,15,250,0.5,200,20e-3]
        if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,20e-3]
    if stat == 'maxe':
        maxes = [60,20,250,0.5,200,2.0]
        if do_rel: maxesx = [0.6,0.02,0.6,0.3,0.6,2.0]
    if stat == 'kgem':
        maxes = [0.9,0.9,0.9,0.9,0.9,20e-3]
    #maxes = [0.0,0.0,0.0,0.0,0.0,0.0] # turn off zoom

    plot_vars = [plot_vars[i] for i in use_vars]
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
    leg_titl = [leg_titl[i] for i in use_vars]
    leg_titl0 = [leg_titl0[i] for i in use_vars]
    leg_titlm = [leg_titlm[i] for i in use_vars]
    maxes = [maxes[i] for i in use_vars]

    method_name = [method_name[i] for i in use_meth]
    plt_name = [plt_name[i] for i in use_meth]
    method_name0 = np.copy(method_name)
    plt_name0 = np.copy(plt_name)

    if do_heat:
        ncol = len(use_meth)
        nrow = len(plot_vars)//ncol
    else:
        ncol = 2
        nrow = len(plot_vars)//ncol + 1

    fig_fil = 'Hrly_diff_scat_{}_{}'
    if do_rel: fig_fil = fig_fil + '_rel'
    if do_heat:
        fig_fil = '{}'+fig_fil + '_heat'
        if sum(maxes)>0: fig_fil = fig_fil + '_zoom'
        fig_fil = fig_fil.format(','.join(method_name),','.join(plot_vars),stat)
    else:
        fig_fil = fig_fil.format(','.join(plot_vars),stat)
    fig_fil = fig_fil + '_compressed.png'

    if 'compressed' in fig_fil: 
        if do_heat: 
            plt.rcParams.update({'font.size': 33})
        else:
            plt.rcParams.update({'font.size': 27})
        fig,axs = plt.subplots(nrow,ncol,figsize=(17*ncol,13*nrow),constrained_layout=do_heat)
    else:
        if do_heat: 
            plt.rcParams.update({'font.size': 120})
        else:
            plt.rcParams.update({'font.size': 100})    
    if not do_heat: fig.subplots_adjust(hspace=0.33, wspace=0.17) # Adjust the bottom margin, vertical space, and horizontal space
        
    for i,(var,lx,ly,plt_t,leg_t,leg_t0,leg_tm,rep,mx) in enumerate(zip(plot_vars,logx,logy,plt_titl,leg_titl,leg_titl0,leg_titlm,rep,maxes)): 
        run_loop(i,var,lx,ly,plt_t,leg_t,leg_t0,leg_tm,rep,mx)

    # Remove the extra subplots
    if not do_heat and (len(plot_vars)) < nrow*ncol:
        for i in range(len(plot_vars),nrow*ncol):
            r = i//ncol
            c = i-r*ncol
            fig.delaxes(axs[r, c])

    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)


if do_balance:
 
 # Specify variables of interest
    if do_heat:
        use_vars = [0,1,2]
        use_meth = [0,1,3]
        logx = np.ones(len(use_vars)) # log scale x axis
        logy = np.ones(len(use_vars)) # log scale y axis
    else:
        use_vars = [0,1,2,3]
        use_meth = [0,1,2,3,4,5]
        logx = np.zeros(len(use_vars)) # no log scale x axis
        logx = np.ones(len(use_vars)) # log scale x axis
        logy = np.ones(len(use_vars)) # log scale y axis

    rep = np.zeros(len(use_vars)) 
    if do_heat:
        use_vars = [val for val in use_vars for _ in range(len(use_meth))]
        logy = [val for val in logy for _ in range(len(use_meth))]
        logx = [val for val in logx for _ in range(len(use_meth))]
        rep = [int(val+_) for val in rep for _ in range(len(use_meth))]

    plot_vars = ['balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceCasNrg','wallClockTime']
    comp_vars = ['balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','numberFluxCalc']
    plt_titl = ['vegetation balance','snow balance','soil balance', 'canopy air space and aquifer balance', 'wall clock time']
    leg_titl = ['$W~m^{-3}$'] * 4 + ['$s$']
    leg_titl0 =['$kg~m^{-2}~s^{-1}$'] * 4 + ['$num$']
    if fixed_Mass_units: leg_titl0 = ['s^{-1}$'] * 3 + ['m~s^{-1}$'] + ['$num$']

    plot_vars = [plot_vars[i] for i in use_vars]
    comp_vars = [comp_vars[i] for i in use_vars]
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
    leg_titl = [leg_titl[i] for i in use_vars]
    leg_titl0 = [leg_titl0[i] for i in use_vars]
    method_name2 = [method_name2[i] for i in use_meth]
    plt_name2 = [plt_name2[i] for i in use_meth]

    method_name20 = np.copy(method_name2)
    plt_name20 = np.copy(plt_name2)

    if do_heat:
        ncol = len(use_meth)
        nrow = len(plot_vars)//ncol
    else:
        ncol = 2
        nrow = len(plot_vars)//ncol + 1

    fig_fil = 'Hrly_balance_scat_{}'
    if do_heat:
        fig_fil = '{}'+fig_fil + '_heat'
        fig_fil = fig_fil.format(','.join(method_name2),','.join(plot_vars),stat)
    else:
        fig_fil = fig_fil.format(','.join(plot_vars),stat)
    fig_fil = fig_fil + '_compressed.png'

    if 'compressed' in fig_fil: 
        if do_heat: 
            plt.rcParams.update({'font.size': 33})
        else:
            plt.rcParams.update({'font.size': 27})
        fig,axs = plt.subplots(nrow,ncol,figsize=(17*ncol,13*nrow),constrained_layout=do_heat)
    else:
        if do_heat: 
            plt.rcParams.update({'font.size': 120})
        else:
            plt.rcParams.update({'font.size': 100})
        fig,axs = plt.subplots(nrow,ncol,figsize=(70*ncol,54*nrow),constrained_layout=do_heat)   
    if not do_heat: fig.subplots_adjust(hspace=0.33, wspace=0.17) # Adjust the bottom margin, vertical space, and horizontal space

    for i,(var,comp,lx,ly,leg_t,leg_t0,plt_t,rep) in enumerate(zip(plot_vars,comp_vars,logx,logy,leg_titl,leg_titl0,plt_titl,rep)): 
        run_loopb(i,var,comp,lx,ly,leg_t,leg_t0,plt_t,rep)

    # Remove the extra subplots
    if not do_heat and (len(plot_vars)) < nrow*ncol:
        for i in range(len(plot_vars),nrow*ncol):
            r = i//ncol
            c = i-r*ncol
            fig.delaxes(axs[r, c])
    # Save
    plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
