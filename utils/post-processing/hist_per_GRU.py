# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python hist_per_GRU.py [stat]
# where stat is rmse or maxe or kgem or rmnz or avge

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy
import pandas as pd

do_box = True # true is plot boxplot instead of CDF/histogram
do_rel = False # true is plot relative to the benchmark simulation
do_hist = False # true is plot histogram instead of CDF
run_local = True # true is run on local machine, false is run on cluster
fix_units_soil = False # true is convert to storage units, only works for soil because of known and constant depth in default setup
fix_wall_actors = False # true then scale reference solution for wall clock time
comp_wall_actors_plot = False # true then plot the wall clock time comparison
comp_wall_event_plot = False # true then plot the event detection time comparison
no_snow = False # true is only plot snow free simulations
# these options are for the boxplot only
showfliers = False # true is show outliers in boxplot
do_violin = False # true is plot violin plot instead of boxplot
vio_points = 10000 # number of points to consider in kernel estimation of violin plot, bigger is better but slow (100 default, 10000 is good for all of N.America)

if run_local: 
    stat = 'avge'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path(os.path.expanduser('~/statistics'))
    

#method_name=['be1','be16','be32','sun6'] #maybe make this an argument
#plt_name=['BE1','BE16','BE32','SUNDIALS'] #maybe make this an argument
#method_name=['sun5cm_noev','sun5cm_ev','sun5en_noev','sun5en_ev','sun8en_noev'] 
#plt_name=['SUNDIALS temp no ev','SUNDIALS temp', 'SUNDIALS enth no ev', 'SUNDIALS enth','reference soln no ev']
method_name=['be8','be8cm','be8en','sun5cm','sun5en'] 
plt_name=['BE8 common','BE8 temp','BE8 mixed','SUNDIALS temp', 'SUNDIALS enth']
#method_name2=method_name
#plt_name2=plt_name
#method_name2=method_name +['sun8en_ev']
#method_name2=method_name +['sun8enOrigWall']
method_name2=method_name +['sun8en']
plt_name2=plt_name +['reference soln']
method_name3=method_name[0:3]
plt_name3=plt_name[0:3]

num_bins = 1000
auto_col = plt.rcParams['axes.prop_cycle'].by_key()['color']

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Define the power transformation function
def power_transform(x):
    return x ** 0.5  # Adjust the exponent as needed

# Simulation statistics file locations
use_vars = []
rep = [] # mark the repeats
#use_vars = [4,4,1,1]
#rep = [1,2,1,2] # mark the repeats
settings0= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','scalarRootZoneTemp']
settings = [settings0[i] for i in use_vars]

use_vars2 = [3,3]
rep2 = [1,2] # mark the repeats
use_vars2 = [8]
rep2 = [0] # mark the repeats
#use_vars2 = [3,3]
#rep2 = [1,2] # mark the repeats
#use_vars2 = [1,1,2,2,3,3]
#rep2 = [1,2,1,2,1,2] # mark the repeats
settings20= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']
settings2 = [settings20[i] for i in use_vars2]

use_vars3 = []
rep3 = [] # mark the repeats
#use_vars3 = [0,1,2,3,0,1,2,3]
#rep3 = [1,1,1,1,2,2,2,2] # mark the repeats
settings30= ['numberStateSplit','numberDomainSplitNrg','numberDomainSplitMass','numberScalarSolutions','meanStepSize']
settings3 = [settings30[i] for i in use_vars3]

viz_fil = method_name.copy()
viz_fl2 = method_name2.copy()
viz_fl3 = method_name3.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_accuracy.nc'
for i, m in enumerate(method_name2):
    viz_fl2[i] = m + '_hrly_diff_bals_balance.nc'
for i, m in enumerate(method_name3):
    viz_fl3[i] = m + '_hrly_diff_steps_split.nc'

# Specify variables of interest
plot_vars = settings.copy()
plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','top 3m soil temperature']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$K$']
if (len(use_vars)>1): 
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
else:
    plt_titl = [f"{plt_titl[i]}" for n,i in enumerate(use_vars)]
leg_titl = [leg_titl[i] for i in use_vars]

plot_vars2 = settings2.copy()
plt_titl2 = ['canopy air space enthalpy balance','vegetation enthalpy balance','snow enthalpy balance','soil enthalpy balance','vegetation mass balance','snow mass balance','soil mass balance','aquifer mass balance', 'wall clock time']
leg_titl2 = ['$W~m^{-3}$'] * 4 + ['$kg~m^{-3}~s^{-1}$'] * 3 + ['$kg~m^{-2}~s^{-1}$']+ ['$s$']
if fix_units_soil: 
    leg_titl2[3] = ['$kJ~m^{-2}$'] 
    leg_titl2[6] = ['$kg~m^{-2}']
if (len(use_vars)+len(use_vars2)>1): 
    plt_titl2 = [f"({chr(97+n + len(use_vars))}) {plt_titl2[i]}" for n,i in enumerate(use_vars2)]
else:
    plt_titl2 = [f"{plt_titl2[i]}" for n,i in enumerate(use_vars2)]
leg_titl2 = [leg_titl2[i] for i in use_vars2]

plot_vars3 = settings3.copy()
plt_titl3 = ['number of state splits','number of energy domain splits','number of mass domain splits','number of scalar solutions','mean step size']
leg_titl3 = [''] * 4 + ['$s$']
if (len(use_vars)+len(use_vars2)+len(use_vars3)>1): 
    plt_titl3 = [f"({chr(97+n + len(use_vars)+len(use_vars2))}) {plt_titl3[i]}" for n,i in enumerate(use_vars3)]
else:
    plt_titl3 = [f"{plt_titl3[i]}" for n,i in enumerate(use_vars3)]
leg_titl3 = [leg_titl3[i] for i in use_vars3]

if do_box:
    fig_fil = 'Hrly_diff_box_{}_{}_zoom'
else:
    if do_hist:
        fig_fil = 'Hrly_diff_hist_{}_{}_zoom'
    else:
        fig_fil = 'Hrly_diff_cdf_{}_{}_zoom'
        #if len(use_vars3)>0: fig_fil = 'Hrly_diff_cdf_{}_{}'
if do_rel: fig_fil = fig_fil+'_rel'
if no_snow: fig_fil = fig_fil + '_nosnow'
fig_fil = fig_fil +'_compressed.png'
if len(use_vars)>0: 
    fig_fil = fig_fil.format('accuracy',stat)
elif len(use_vars2)>0: # and len(use_vars)==0:
    fig_fil = fig_fil.format('wallclock','mean')
elif len(use_vars3)>0: 
    fig_fil = fig_fil.format('split','mean')

maxes_m = [99,15,99,99,6]
if do_rel: maxes_m = [0.4,0.007,0.6,0.15,0.0015]
if stat == 'avge':
    stat2 = 'mean'
    maxes = [99,7,99,99,0.28]
    if do_rel: maxes = [0.4,0.007,0.6,0.15,0.0015]
if stat == 'rmse' or stat=='rmnz':
    stat2 = 'mean'
    maxes = [2,15,250,0.08,200]
    if do_rel: maxes = [0.4,0.007,0.6,0.15,0.0015]
if stat == 'maxe':
    stat2 = 'amax'
    if stat == 'maxe': maxes = maxes_m
if stat == 'kgem':
    stat2 = 'mean'
    maxes = [0.9,0.9,0.9,0.9,0.9]
maxes = [maxes[i] for i in use_vars]
for i in range(len(maxes)):
    #if rep[i]==2: maxes[i] = maxes[i]*2.5 #clunky way to increase the plot_range for the second repeat
    if rep[i]==2: maxes[i] = maxes_m[use_vars[i]] #clunky way to increase the plot_range for the second repeat

if stat2 == 'mean':
    maxes2 = [5e1,5e1,5e1,5e1]+[1e-7,1e-5,1e-7,1e-8] + [2e-2]
if stat2 == 'amax':
    maxes2 = [5e3,5e3,5e3,5e3]+[1e-5,1e-3,1e-5,1e-6] + [2.0]
maxes2 = [maxes2[i] for i in use_vars2]
for i in range(len(maxes2)):
    if rep2[i]==2: maxes2[i] = maxes2[i]*1e2 #clunky way to increase the plot_range for the second repeat

stat3 = 'mean'
maxes3 = [1e2,1e2,1e2,1e2,1e-7]
maxes3 = [maxes3[i] for i in use_vars3]

summa = {}
summa1 = {}
summa2 = {}
if len(use_vars)>0:
    for i, m in enumerate(method_name):
        # Get the aggregated statistics of SUMMA simulations
        summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
if len(use_vars2)>0:
    for i, m in enumerate(method_name2):
        summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])
    if (fix_wall_actors or comp_wall_actors_plot) and 'wallClockTime' in settings2:
        summa1['be8Old'] = xr.open_dataset(viz_dir/'be8NrgOld_hrly_diff_bals_balance.nc')
        summa1['sun5enOld'] = xr.open_dataset(viz_dir/'sun5enNrgOld_hrly_diff_bals_balance.nc')

if len(use_vars3)>0:
    for i, m in enumerate(method_name3):
        summa2[m] = xr.open_dataset(viz_dir/viz_fl3[i])

if no_snow:
    summa[method_name[0]] = xr.open_dataset(viz_dir/viz_fil[0]) # will be a problem if this does not exist
    if len(use_vars)>0:
        for m in method_name:
            summa[m] = summa[m].where(summa[method_name[0]]['scalarSWE'].sel(stat='mean_ben') == 0)
    if len(use_vars2)>0:
        for m in method_name2:
            summa1[m] = summa1[m].where(summa[method_name[0]]['scalarSWE'].sel(stat='mean_ben') == 0)
        
    if len(use_vars3)>0:
        for m in method_name3:
            summa2[m] = summa2[m].where(summa[method_name[0]]['scalarSWE'].sel(stat='mean_ben') == 0)

    
##Figure

plt.rcParams['xtick.color'] = 'black'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['ytick.major.width'] = 2
# fix size for now
ncol = 2
nrow = 3

if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 27})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(nrow,ncol,figsize=(17*ncol,17*nrow))
else:
    fig,axs = plt.subplots(nrow,ncol,figsize=(70*ncol,80*nrow))
fig.subplots_adjust(hspace=0.2, wspace=0.12) # Adjust the bottom margin, vertical space, and horizontal space
#fig.suptitle('Histograms of Hourly Statistics for each GRU', fontsize=40,y=1.0)
    
def run_loop(i,var,mx,rep,stat):
    r = i//ncol
    c = i-r*ncol
    if rep == 1: stat = 'avge'
    if rep == 2: stat = 'maxe'
    stat0 = stat
    if stat == 'rmse' or stat == 'kgem' or stat == 'avge': 
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
        s_rel = summa[method_name[0]][var].sel(stat=statr)
        for m in method_name:
            s = summa[m][var].sel(stat=stat0)
            if do_rel and var != 'wallClockTime': s = s/s_rel
            if stat == 'maxe': s = np.fabs(s) # make absolute value norm
            mx = max(s.max(),mx)
            if stat == 'kgem': mn = min(s.min(),mn)

    # Data
    s_rel = summa[method_name[0]][var].sel(stat=statr)
    for m in method_name:
        s = summa[m][var].sel(stat=stat0)
        if do_rel and var != 'wallClockTime': s = s/s_rel  
        if stat == 'maxe': s = np.fabs(s) # make absolute value norm
        plot_range = (0,mx)
        if stat=='kgem' and var!='wallClockTime': 
            plot_range = (mn,1)
        elif var=='wallClockTime':
            plot_range = (0.0008,mx)
        if do_box:
            data = np.fabs(s.values)
            data = data[~np.isnan(data)]
            if do_violin:
                vplot = axs[r, c].violinplot(dataset=[data],positions=[len(method_name) - method_name.index(m)],vert=False,showextrema=showfliers,points=vio_points)
                for pc in vplot['bodies']:
                    pc.set_facecolor(auto_col[method_name.index(m)])
                    pc.set_edgecolor('black')
                    pc.set_alpha(1)
            else:
                axs[r, c].boxplot(data,vert=False, positions=[len(method_name) - method_name.index(m)], widths=0.6,patch_artist=True,medianprops=dict(color='black'),boxprops=dict(facecolor=auto_col[method_name.index(m)]),showfliers=showfliers)

        else:
            if do_hist: 
                np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=3.0,range=plot_range)
            else: #cdf
                sorted_data = np.sort(np.fabs(s))
                valid_data = sorted_data[~np.isnan(sorted_data)]
                yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
                axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=3.0)
        print("max, min, mean without nans", s.where(lambda x: ~np.isnan(x)).max().values, s.where(lambda x: ~np.isnan(x)).min().values, s.where(lambda x: ~np.isnan(x)).mean().values, m, var)

    if stat0 == 'rmse': stat_word = 'RMSE'
    if stat0 == 'rmnz': stat_word = 'RMSE' # no 0s'
    if stat0 == 'maxe': stat_word = 'max abs error'
    if stat0 == 'kgem': stat_word = 'KGE"'
    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'mnnz': stat_word = 'mean' # no 0s'
    if stat0 == 'amax': stat_word = 'max'
    if stat0 == 'avge': stat_word = 'mean abs error'
    
    if statr == 'mean_ben': statr_word = 'mean'
    if statr == 'mnnz_ben': statr_word = 'mean' # no 0s'
    if statr == 'amax_ben': statr_word = 'max'
    
    if c==0 and not do_box: axs[r,c].legend(plt_name)
    titl = plt_titl[i]
    if no_snow: titl = titl + ' (snow-free GRUs)'
    if rep>0: titl = titl #+ ' '+ stat_word
    axs[r,c].set_title(titl)
    if stat=='rmse' or stat=='rmnz' or stat=='maxe' or stat=='mean' or stat=='avge': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl[i]))
    if stat=='kgem': axs[r,c].set_xlabel(stat_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    if do_box:
        axs[r,c].set_xlim(plot_range)
        axs[r, c].set_ylabel('')
        axs[r, c].set_yticklabels('')
        axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
        if mx<1: # Rotate x-axis labels
             axs[r, c].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot
        axs[r, c].set_yticks(range(1, len(method_name) + 1))
        if(c==0): axs[r, c].set_yticklabels(plt_name[::-1])
    else:
        if do_hist: 
            axs[r,c].set_ylabel('GRU count')
            if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
        else:
            axs[r,c].set_xlim(plot_range)
            axs[r,c].set_ylabel('cumulative distribution')
            if(c>=1): axs[r, c].set_ylabel('')
            axs[r,c].set_ylim([0.0, 1.0])
            axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
            if mx<1: # Rotate x-axis labels
                axs[r,c].tick_params(axis='x', rotation=45)


def run_loopb(i,var,mx,rep,stat2):
    r = (i+len(use_vars))//ncol
    c = (i+len(use_vars))-r*ncol
    stat0 = np.copy(stat2)
    if rep == 1: stat0 = 'mean'
    if rep == 2: stat0 = 'amax'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx*1e-9
        if any(substring in var for substring in ['VegNrg', 'SnowNrg', 'SoilNrg']):
            mn = mx*1e-9
        if var=='wallClockTime': mn = 0.0008
        if fix_units_soil and 'Soil' in var:
            mn = mn*3600*4.0 # mult by time step and depth to get storage
            mx = mx*3600*4.0
            if 'Nrg' in var:
                mn=mn*1e-3
                mx=mx*1e-3
    else:
        mx = 0.0
        mn = 1.0
        for m in method_name2:
            # Get the statistics, remove 9999 (should be nan, but just in case)
            s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
            if var=='wallClockTime': s = s.where(lambda x: x != 0) # Actors simulations may have 0
            mx = max(s.max(),mx)
            mn = min(s.min(),mn)
    # Data
    combined_s2 = []
    combined_s_saved = []
    for m in method_name2:
        s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)            
        if var=='wallClockTime': s = s.where(lambda x: x != 0) # water bodies should be 0
        if fix_units_soil and 'Soil' in var: 
            s = s*3600*4.0 # mult by time step and depth to get storage
            if 'Nrg' in var: s = s*1e-3

        plot_range = (mn,mx)

        if (fix_wall_actors or comp_wall_actors_plot) and 'wallClockTime' in var: 
            from scipy.stats import linregress
            if m in ['be8', 'sun5en']:
                s_saved = s
                s2 = summa1[f'{m}Old'][var].sel(stat=stat0).where(lambda x: x != 9999)
                s2 = s2.where(lambda x: x != 0)  # water bodies should be 0
                mask = ~np.isnan(s2.values) & ~np.isnan(s_saved.values)
                s2 = s2[mask]
                s_saved = s_saved[mask]
                combined_s2.append(s2.values)
                combined_s_saved.append(s_saved.values)
                first_len = len(s2)

            if m=='sun8en': # assumes sun8en is the last one
                combined_s2 = np.concatenate(combined_s2)
                combined_s_saved = np.concatenate(combined_s_saved)
                # Least squares fit
                A = combined_s2[:, np.newaxis]
                fac, _, _, _ = np.linalg.lstsq(A, combined_s_saved, rcond=None)
                fac = fac[0]
                print(f'Best fit least squares ratio (slope={fac:.4f})')
                slope, intercept, r_value, p_value, std_err = linregress(combined_s2, combined_s_saved)
                print(f'Best fit regression line (slope={slope:.4f}, intercept={intercept:.4f}, corr coeff={r_value:.2e})')
                # slope = 0.8731
                # intercept = -0.0003
                # corr coeff = 0.880
                # x is s2, non-actors Graham way (Anvil is faster)
                s = s * slope + intercept
                
        if comp_wall_event_plot and 'wallClockTime' in var:
            from scipy.stats import linregress
            if m in ['sun5cm_ev','sun5en_ev','sun8en_ev']:
                s_saved = s
                s2 = summa1[f'{m[:-3]}_noev'][var].sel(stat=stat0).where(lambda x: x != 9999)
                s2 = s2.where(lambda x: x != 0)  # water bodies should be 0
                mask = ~np.isnan(s2.values) & ~np.isnan(s_saved.values)
                s2 = s2[mask]
                s_saved = s_saved[mask]
                combined_s2.append(s2.values)
                combined_s_saved.append(s_saved.values)
                first_len = len(s2)

            if m=='sun8en_ev': # assumes sun8en is the last one
                combined_s2 = np.concatenate(combined_s2)
                combined_s_saved = np.concatenate(combined_s_saved)
                # Least squares fit
                A = combined_s2[:, np.newaxis]
                fac, _, _, _ = np.linalg.lstsq(A, combined_s_saved, rcond=None)
                fac = fac[0]
                print(f'Best fit least squares ratio (slope={fac:.4f})')
                slope, intercept, r_value, p_value, std_err = linregress(combined_s2, combined_s_saved)
                print(f'Best fit regression line (slope={slope:.4f}, intercept={intercept:.4f}, corr coeff={r_value:.2e})')
                # slope = 1.0423
                # intercept = -0.0001
                # corr coeff = 0.987
                # x is s2, no_ev way (with event detection is slower)
                # note, do not apply the fit to the event detection time

        if do_box:
            data = np.fabs(s.values)
            data = data[~np.isnan(data)]
            if do_violin:
                vplot = axs[r, c].violinplot(dataset=[data],positions=[len(method_name2) - method_name2.index(m)],vert=False,showextrema=showfliers,points=vio_points)
                for pc in vplot['bodies']:
                    pc.set_facecolor(auto_col[method_name2.index(m)])
                    pc.set_edgecolor('black')
                    pc.set_alpha(1)
            else:
                axs[r, c].boxplot(data,vert=False, positions=[len(method_name2) - method_name2.index(m)], widths=0.6,patch_artist=True,medianprops=dict(color='black'),boxprops=dict(facecolor=auto_col[method_name2.index(m)]),showfliers=showfliers)
        else:
            if do_hist: 
                np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=3.0,range=plot_range)
            else: #cdf
                sorted_data = np.sort(np.fabs(s))
                valid_data = sorted_data[~np.isnan(sorted_data)]
                yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
                axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=3.0)
        print("max, min, mean without nans", s.where(lambda x: ~np.isnan(x)).max().values, s.where(lambda x: ~np.isnan(x)).min().values, s.where(lambda x: ~np.isnan(x)).mean().values, m, var)

    if stat0 == 'mean': 
        if var == 'wallClockTime': 
            stat_word = 'mean'
        else:
            stat_word = 'mean abs balance'
    if stat0 == 'amax': 
        if var == 'wallClockTime': 
            stat_word = 'max'
        else:
            stat_word = 'max abs balance'

    if c==0 and not do_box: axs[r,c].legend(plt_name2)
    titl = plt_titl2[i]
    if no_snow: titl = titl + ' (snow-free GRUs)'
    if rep>0: titl = titl #+ ' '+ stat_word
    axs[r,c].set_title(titl)
    axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl2[i]))   

    if do_box:
        axs[r,c].set_xlim(plot_range)
        axs[r, c].set_ylabel('')
        axs[r, c].set_yticklabels('')
        axs[r,c].set_xscale('log') #log x axis
        if var=='wallClockTime': 
            axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
            axs[r, c].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot
        axs[r, c].set_yticks(range(1, len(method_name2) + 1))
        if(c==0): axs[r, c].set_yticklabels(plt_name2[::-1])
    else:
        if do_hist: 
            axs[r,c].set_ylabel('GRU count')
            if(c==1): axs[r, c].set_ylabel('')
            if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
        else:
            axs[r,c].set_xlim(plot_range)
            axs[r,c].set_ylabel('cumulative distribution')
            if(c>=1): axs[r, c].set_ylabel('')
            axs[r,c].set_ylim([0.0, 1.0])
            axs[r,c].set_xscale('log') #log x axis
            if var=='wallClockTime': 
                axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
                axs[r, c].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot
    
    if comp_wall_actors_plot:
        fig.subplots_adjust(hspace=0.2, wspace=0.2) # Adjust the bottom margin, vertical space, and horizontal space
        axs[r, c + 1].scatter(combined_s2[first_len:], combined_s_saved[first_len:], alpha=0.5, color=auto_col[4], label='SUNDIALS enth')
        axs[r, c + 1].scatter(combined_s2[:first_len], combined_s_saved[:first_len], alpha=0.5, color=auto_col[0], label='BE8 common')
        axs[r, c+1].set_xlabel('Graham time [s]')
        axs[r, c+1].set_ylabel('Anvil Actors time [s]')
        axs[r, c+1].set_title('wall clock time comparison')
        axs[r, c+1].set_xlim(combined_s_saved.min(),combined_s2.max()) 
        axs[r, c+1].set_ylim(combined_s_saved.min(),combined_s2.max())
        axs[r, c+1].plot(combined_s2, intercept + slope * combined_s2, color='black',linewidth=3.0)
        axs[r, c+1].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot

    if comp_wall_event_plot:
        fig.subplots_adjust(hspace=0.2, wspace=0.2) # Adjust the bottom margin, vertical space, and horizontal space
        axs[r, c + 1].scatter(combined_s2[first_len:], combined_s_saved[first_len:], alpha=0.5, color=auto_col[5], label='reference soln')
        axs[r, c + 1].scatter(combined_s2[:first_len], combined_s_saved[:first_len], alpha=0.5, color=auto_col[3], label='SUNDIALS enth')
        axs[r, c+1].set_xlabel('no event detection time [s]')
        axs[r, c+1].set_ylabel('event detection time time [s]')
        axs[r, c+1].set_title('wall clock time comparison')
        axs[r, c+1].set_xlim(combined_s_saved.min(),combined_s2.max()) 
        axs[r, c+1].set_ylim(combined_s_saved.min(),combined_s2.max())
        axs[r, c+1].plot(combined_s2, intercept + slope * combined_s2, color='black',linewidth=3.0)
        axs[r, c+1].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot


def run_loop3(i,var,mx,rep,stat3):
    r = (i+len(use_vars)+len(use_vars2))//ncol
    c = (i+len(use_vars)+len(use_vars2))-r*ncol
    stat0 = np.copy(stat3)
    if rep == 1: stat0 = 'mean'
    if rep == 2: stat0 = 'amax'

    mx = 0.0
    mn = 1.0
    for m in method_name3:
        s = summa2[m][var].sel(stat=stat0)
        mx = max(s.max(),mx)
        mn = min(s.min(),mn)

    # Data
    for m in method_name3:
        s = summa2[m][var].sel(stat=stat0)
        plot_range = (0,mx)
        if do_box:
            data = np.fabs(s.values)
            data = data[~np.isnan(data)]
            if do_violin:
                vplot = axs[r, c].violinplot(dataset=[data],positions=[len(method_name3) - method_name3.index(m)],vert=False,showextrema=showfliers,points=vio_points)
                for pc in vplot['bodies']:
                    pc.set_facecolor(auto_col[method_name3.index(m)])
                    pc.set_edgecolor('black')
                    pc.set_alpha(1)
            else:
                axs[r, c].boxplot(data,vert=False, positions=[len(method_name3) - method_name3.index(m)], widths=0.6,patch_artist=True,medianprops=dict(color='black'),boxprops=dict(facecolor=auto_col[method_name3.index(m)]),showfliers=showfliers)
        else:
            if do_hist: 
                np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=3.0,range=plot_range)
            else: #cdf
                sorted_data = np.sort(np.fabs(s))
                valid_data = sorted_data[~np.isnan(sorted_data)]
                yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
                axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=3.0)
        print("max, min, mean without nans", s.where(lambda x: ~np.isnan(x)).max().values, s.where(lambda x: ~np.isnan(x)).min().values, s.where(lambda x: ~np.isnan(x)).mean().values, m, var)

    if stat0 == 'mean': stat_word = 'mean per data window'
    if stat0 == 'amax': stat_word = 'max per data window'

    if c==0 and not do_box: axs[r,c].legend(plt_name3)
    titl = plt_titl3[i]
    if no_snow: titl = titl + ' (snow-free GRUs)'
    if rep>0: titl = titl #+ ' '+ stat_word
    axs[r,c].set_title(titl)
    axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl3[i]))   

    if do_box:
        axs[r,c].set_xlim(plot_range)
        axs[r, c].set_ylabel('')
        axs[r, c].set_yticklabels('')
        #axs[r,c].set_xscale('log') #log x axis
        axs[r, c].set_yticks(range(1, len(method_name3) + 1))
        if(c==0): axs[r, c].set_yticklabels(plt_name3[::-1])
    else:
        if do_hist: 
            axs[r,c].set_ylabel('GRU count')
            if(c==1): axs[r, c].set_ylabel('')
        else:
            axs[r,c].set_xlim(plot_range)
            axs[r,c].set_ylabel('cumulative distribution')
            if(c>=1): axs[r, c].set_ylabel('')
            axs[r,c].set_ylim([0.0, 1.0])
            #axs[r,c].set_xscale('log') #log x axis
        if 'zoom' in fig_fil:
            axs[r,c].set_ylim([0.98, 1.0])
 

if len(use_vars) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars,maxes,rep)): 
        run_loop(i,var,mx,rep,stat)
if len(use_vars2) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars2,maxes2,rep2)): 
        run_loopb(i,var,mx,rep,stat2)
if len(use_vars3) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars3,maxes3,rep3)): 
        run_loop3(i,var,mx,rep,stat3)


# Remove the extra subplots
if (len(plot_vars)+len(plot_vars2)+len(plot_vars3)) < ncol*nrow:
    for i in range((len(plot_vars)+len(plot_vars2)+len(plot_vars3)),ncol*nrow):
        r = i//ncol
        c = i-r*ncol
        if (r==0 and c==1 and comp_wall_actors_plot): continue
        if (r==0 and c==1 and comp_wall_event_plot): continue
        fig.delaxes(axs[r, c])

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
