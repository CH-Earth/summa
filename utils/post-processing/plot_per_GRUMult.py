# written originally by W. Knoben, modified by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    Catchment shapefile with GRU delineation
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# To improve visualization of large lakes, HydroLAKES lake delineations are plotted on top of the catchment GRUs and river segments.
# Dealing with HydroLAKES inputs is not considered within scope of the workflow and therefore requires some manual downloading and preprocessing of this data for those who wish to reproduce this step.
# The relevant code is easily disabled by switching the plot_lakes = True flag to False.

# Run:
# python plot_per_GRUMult.py [stat]
# where stat is rmse or maxe or kgem or mean or amax

# modules
import sys
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy
import pyproj
import fiona
import geopandas as gpd
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.ticker import ScalarFormatter


do_rel = True # true is plot relative to the benchmark simulation
one_plot = False # true is one plot, false is multiple plots (one per variable)
run_local = False # true is run on local machine (only does testing), false is run on cluster
more_mean = False # true is plot mean/amax extra variables in a balance file

if run_local: 
    stat = 'mean'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')


# NOTE: method_name 'ref' will plot the reference solution, 'diff' will plot the difference between two simulations
#       method_name 'diff' requires the specification of the two simulations to subtract in the variables from_meth and sub_meth

method_name=['be1','be16','be32','sundials_1en6','ref']
plt_name0=['SUMMA-BE1','SUMMA-BE16','SUMMA-BE32','SUMMA-SUNDIALS','reference solution']
plt_nameshort=plt_name0
method_name=['be1','be1cm','be1en','sundials_1en6cm','sundials_1en6en','diff','ref']
plt_name0=['BE1 common thermo. eq.','SUMMA-BE1 temperature thermo. eq.','SUMMA-BE1 mixed thermo. eq.','SUMMA-SUNDIALS temperature thermo. eq.','SUMMA-SUNDIALS enthalpy thermo. eq.','SUMMA-BE1 common - mixed','reference solution']
plt_nameshort=['BE1 common','BE1 temp','BE1 mixed','SUNDIALS temp','SUNDIALS enth','BE1 common - mixed','reference soln']

if one_plot: plt_name0 = plt_nameshort

from_meth = 'be1' # name of the first simulation in the difference simulation, only used if a method_name is 'diff'
sub_meth = 'be1en' # name of the simulation to subtract in the difference simulation, only used if a method_name is 'diff'

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']

viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
nbatch_hrus = 518 # number of HRUs per batch
if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

if more_mean: # extra vars in a balance file
    plt_titl_exVar = ['rain plus melt','top 4m soil temperature','air temperature','snow water equivalent']
    plot_vars_exVar = ['scalarRainPlusMelt','scalarRootZoneTemp','airtemp','scalarSWE']
    viz_file_exVar = 'exVar_hrly_diff_bals_balance.nc'
    plt_name0_exVar = 'SUMMA-BE1 temperature thermo. eq.'
    plt_nameshort_exVar = 'BE1 temp' # identify method here
    leg_titl_exVar = ['$mm~y^{-1}$','$K$','$K$','$kg~m^{-2}$']
    maxes_exVar = [3000,290,290,100]
    if one_plot: plt_name0_exVar = plt_nameshort_exVar

# Specify variables in files
plot_vars = settings.copy() + ['scalarSWE']
plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','average routed runoff','wall clock time', 'melt with seasonal snow']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$s$','$kg~m^{-2}$']
calc = [0,0,0,0,0,0,1] # 1 if variable needs to be calculated from other variables
melt_thresh = 1/(0.75) # threshold for melt water calculation (divisor is percentage of year no snow, if only melts once)

fig_fil= '_hrly_diff_stats_{}_compressed.png'
if do_rel: fig_fil = '_hrly_diff_stats_{}_rel_compressed.png'

if stat == 'rmse' or stat == 'rmnz': 
    maxes = [2,15,250,0.08,200,10e-3,2] 
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,10e-3,0.6]
if stat == 'maxe': 
    maxes = [15,25,0.8,2,0.3,0.2,15] #[15,25,25e-5,2,1e-7,0.2]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,0.2,0.6]
if stat == 'kgem': 
    maxes = [0.9,0.9,0.9,0.9,0.9,10e-3,0.9]
if stat == 'mean' or stat == 'mnnz': 
    maxes = [100,1700,2000,8,4000,10e-3,100] #[80,1500,5e-5,8,1e-7,10e-3]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,10e-3,1.1]
if stat == 'amax': 
    maxes = [240,1800,3.5,25,7.5,0.2,240] #[240,1800,1e-3,25,2e-6,0.2]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,0.2,1.1]

# Get simulation statistics
summa = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    if m!='diff' and m!='ref': summa[m] = xr.open_dataset(viz_dir/viz_fil[i])

if more_mean: summa['exVar'] = xr.open_dataset(viz_dir/viz_file_exVar)



# Function to extract a given setting from the control file
def read_from_control( file, setting ):

    # Open controlFile and ...
    with open(file) as contents:
        for line in contents:

            # ... find the line with the requested setting
            if setting in line and not line.startswith('#'):
                break

    # Extract the setting's value
    substring = line.split('|',1)[1]      # Remove the setting's name (split into 2 based on '|', keep only 2nd part)
    substring = substring.split('#',1)[0] # Remove comments, does nothing if no '#' is found
    substring = substring.strip()         # Remove leading and trailing whitespace, tabs, newlines

    # Return this value
    return substring

# Function to specify a default path
def make_default_path(suffix):

    # Get the root path
    rootPath = Path( read_from_control(controlFile,'root_path') )

    # Get the domain folder
    domainName = read_from_control(controlFile,'domain_name')
    domainFolder = 'domain_' + domainName

    # Specify the forcing path
    defaultPath = rootPath / domainFolder / suffix

    return defaultPath

if run_local:
    # Make stubs to check if the plots set up properly
    plot_lakes = False
    plot_rivers = False

    # Create a mock DataFrame
    from shapely.geometry import Point

    s = summa[method_name[0]][plot_vars[0]].sel(stat=stat)
    mock_data = {
        'hm_hruid': np.concatenate(([81029662], s.hru.values[-100:])), #s.hru.values[-100:],  # Example HRU IDs
        'geometry': [Point(x, y) for x, y in zip(range(101), range(101))]  # Simple geometries
    }
    bas_albers = gpd.GeoDataFrame(mock_data, geometry='geometry')    
    hm_hruid = 'hm_hruid'  # Correctly define the variable name in the shapefile
    xmin, ymin, xmax, ymax = bas_albers.total_bounds

else:
    # Get the albers shapes
    main = Path('/home/avanb/projects/rpp-kshook/wknoben/CWARHM_data/domain_NorthAmerica/shapefiles/albers_projection')
    plot_lakes = True
    plot_rivers = False

    # Control file handling
    controlFile = 'plot_control_NorthAmerica.txt'

    # HM catchment shapefile path & name
    hm_catchment_path = read_from_control(controlFile,'catchment_shp_path')
    hm_catchment_name = read_from_control(controlFile,'catchment_shp_name')
    # Specify default path if needed
    if hm_catchment_path == 'default':
        hm_catchment_path = make_default_path('shapefiles/catchment') # outputs a Path()
    else:
        hm_catchment_path = Path(hm_catchment_path) # make sure a user-specified path is a Path()

    # Find the GRU and HRU identifiers
    hm_hruid = read_from_control(controlFile,'catchment_shp_hruid')

    ## River network shapefile location and variable names
    river_network_path = read_from_control(controlFile,'river_network_shp_path')
    river_network_name = read_from_control(controlFile,'river_network_shp_name')
    # Specify default path if needed
    if river_network_path == 'default':
        river_network_path = make_default_path('shapefiles/river_network') # outputs a Path()
    else:
        river_network_path = Path(river_network_path) # make sure a user-specified path is a Path()

    # Find the segment ID
    seg_id = read_from_control(controlFile,'river_network_shp_segid')

    ## Load all shapefiles and project to Albers Conformal Conic and reproject
    acc = 'ESRI:102008' # Set the target CRS

    bas_albers = gpd.read_file(main/'basin.shp')
    xmin, ymin, xmax, ymax = bas_albers.total_bounds

    if plot_rivers: riv_albers = gpd.read_file(main/'river.shp')
    if plot_lakes: lak_albers = gpd.read_file(main/'lakes.shp')



# Match the accummulated values to the correct HRU IDs in the shapefile
hru_ids_shp = bas_albers[hm_hruid].astype(int) # hru order in shapefile
for i,plot_var in enumerate(plot_vars):
    stat0 = stat
    if stat == 'rmse' or stat == 'kgem' or stat == 'mean': 
        if plot_var == 'wallClockTime': stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz' or stat == 'mnnz':
        if plot_var == 'wallClockTime': stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe' or stat == 'amax': 
        if plot_var == 'wallClockTime': stat0 = 'amax'
        statr = 'amax_ben'

    if do_rel: s_rel = np.fabs(summa[method_name[0]][plot_var].sel(stat=statr))
    
    if calc[i]:
        if do_rel: s_rel = s_rel.where(summa[method_name[0]][plot_var].sel(stat='mnnz_ben') > melt_thresh*summa[method_name[0]][plot_var].sel(stat='mean_ben'))


    for m in method_name:
        if m=='diff': 
            s = summa[from_meth][plot_var].sel(stat=stat0) - summa[sub_meth][plot_var].sel(stat=stat0)
        elif m=='ref':
            s = np.fabs(summa[method_name[0]][plot_var].sel(stat=statr))
        else:
            s = np.fabs(summa[m][plot_var].sel(stat=stat0))
        if calc[i]:
            if m=='diff': 
                s_from = summa[from_meth][plot_var].sel(stat=stat0)
                s_from = s_from.where(summa[from_meth][plot_var].sel(stat='mnnz') > melt_thresh*summa[from_meth][plot_var].sel(stat='mean'))
                s_sub  = summa[sub_meth][plot_var].sel(stat=stat0)
                s_sub  = s_sub.where(summa[sub_meth][plot_var].sel(stat='mnnz') > melt_thresh*summa[sub_meth][plot_var].sel(stat='mean'))
                s = s_from - s_sub
            elif m=='ref':
                s =s.where(summa[method_name[0]][plot_var].sel(stat='mnnz_ben') > melt_thresh*summa[method_name[0]][plot_var].sel(stat='mean_ben'))
            else:
                s = s.where(summa[m][plot_var].sel(stat='mnnz') > melt_thresh*summa[m][plot_var].sel(stat='mean'))
        if do_rel and plot_var != 'wallClockTime': s = s/s_rel

        # Replace inf and 9999 values with NaN in the s DataArray
        s = s.where(~np.isinf(s), np.nan).where(lambda x: x != 9999, np.nan)

        if plot_var == 'scalarTotalET' and not do_rel:
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz' or stat=='mean': s = s*31557600 # make annual total
            if stat =='maxe' or stat=='amax': s = s*3600 # make hourly max
        if plot_var == 'averageRoutedRunoff' and not do_rel:
            if stat =='rmse' or stat =='rmnz' or stat=='mnnz' or stat=='mean': s = s*31557600*1000 # make annual total
            if stat =='maxe' or stat=='amax': s = s*3600*1000 # make hourly max

        # Create a new column in the shapefile for each method, and fill it with the statistics
        if calc[i]: 
            plot_var1 = plot_var + '_calc'
        else:
            plot_var1 = plot_var
        bas_albers[plot_var1+m] = np.nan
        hru_ind = [i for i, hru_id in enumerate(hru_ids_shp.values) if hru_id in s.hru.values] # if some missing
        bas_albers.loc[hru_ind, plot_var1+m] = s.sel(hru=hru_ids_shp.values[hru_ind]).values

if more_mean: # extra mean/amax variables
    for i,plot_var in enumerate(plot_vars_exVar):
        stat0 = stat
    
        if stat != 'mean' and stat != 'amax': 
            print('Only mean and amax are supported for extra variables')
            sys.exit()

        m = 'exVar'
        s = np.fabs(summa[m][plot_var].sel(stat=stat0))

        # Replace inf and 9999 values with NaN in the s DataArray
        s = s.where(~np.isinf(s), np.nan).where(lambda x: x != 9999, np.nan)

        if plot_var == 'scalarRainPlusMelt':
            if stat=='mean': s = s*31557600*1000 # make annual total
            if stat=='amax': s = s*3600*1000 # make hourly max

        # Create a new column in the shapefile for each method, and fill it with the statistics
        plot_var1 = plot_var
        bas_albers[plot_var1+m] = np.nan
        hru_ind = [i for i, hru_id in enumerate(hru_ids_shp.values) if hru_id in s.hru.values] # if some missing
        bas_albers.loc[hru_ind, plot_var1+m] = s.sel(hru=hru_ids_shp.values[hru_ind]).values



# Select lakes of a certain size for plotting
if plot_lakes:
    minSize = 1000 # km2
    in_domain = (lak_albers['Country'] == 'Canada') | \
                (lak_albers['Country'] == 'United States of America') | \
                (lak_albers['Country'] == 'Mexico')
    out_domain = (lak_albers['Pour_long'] > -80) & (lak_albers['Pour_lat'] > 65) # Exclude Baffin Island
    large_lakes_albers = lak_albers.loc[(lak_albers['Lake_area'] > minSize) & in_domain & (~out_domain) ]
    lake_col = (8/255,81/255,156/255)
    


# Figure
def run_loop(j,var,the_max):
    stat0 = stat
    if stat == 'rmse' or stat == 'kgem' or stat == 'mean': 
        if var == 'wallClockTime': stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz' or stat == 'mnnz':
        if var == 'wallClockTime': stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe' or stat == 'amax': 
        if var == 'wallClockTime': stat0 = 'amax'
        statr = 'amax_ben'

        
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

    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    my_cmap.set_bad(color='white') #nan color white

    if var!='exVar':
        vmin,vmax = 0, the_max
        if (stat =='mean' or stat=='mnnz') and var=='scalarTotalSoilWat' and not do_rel: vmin,vmax = 600, the_max
        if stat =='amax' and var=='scalarTotalSoilWat' and not do_rel: vmin,vmax = 1000, the_max
        if (stat == 'mean' or stat == 'mnnz' or stat == 'amax') and var!='wallClockTime' and do_rel: vmin,vmax = 0.9, the_max
        norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=0.5)

        if stat =='kgem' and var!='wallClockTime':
            my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno')) # copy the default cmap
            my_cmap.set_bad(color='white') #nan color white
            vmin,vmax = the_max, 1.0
            norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=1.5)

        my_cmap2 = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
        my_cmap2.set_bad(color='white') #nan color white
        vmin,vmax =  -the_max/250,the_max/250,
        norm2 = matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        #norm2 = matplotlib.colors.SymLogNorm(vmin=vmin,vmax=vmax,linthresh=0.01,base=1.1)

        for i,m in enumerate(method_name):
            r = i//ncol + base_row
            c = i - (r-base_row)*ncol

            # Plot the data with the full extent of the bas_albers shape
            if m=='diff':
                bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap2, norm=norm2,zorder=0)
                stat_word0 = stat_word+' difference'
                stat_word2 = stat_word
                plt_nm = plt_name[i]
            elif m=='ref' and var != 'wallClockTime':
                # only plot wallClockTime for the reference solution
                plt_nm =''
            else:
                bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)
                stat_word0 = stat_word
                plt_nm = plt_name[i]
            print(f"{'all HRU mean for '}{var+m:<35}{np.nanmean(bas_albers[var+m].values):<10.5f}{' max: '}{np.nanmax(bas_albers[var+m].values):<10.5f}")
            axs[r,c].set_title(plt_nm)
            axs[r,c].axis('off')
            axs[r,c].set_xlim(xmin, xmax)
            axs[r,c].set_ylim(ymin, ymax)

            # Custom colorbar
            if i==len(method_name)-1:
                if m=='diff':
                    sm = matplotlib.cm.ScalarMappable(cmap=my_cmap2, norm=norm2)
                    sm2 = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
                else:
                    sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
                sm.set_array([])
                if one_plot:
                    if m=='diff': # only works if diff is last on list
                        cbr = fig.colorbar(sm, ax=axs_list[r*len(method_name):(r+1)*len(method_name)],aspect=27/nrow,location='right')
                        cbr2 = fig.colorbar(sm2, ax=axs_list[(r+1)*len(method_name)-1:(r+1)*len(method_name)],aspect=27/nrow,location='left')
                        cbr2.ax.yaxis.set_ticks_position('right')
                        cbr2.ax.yaxis.set_label_position('right')
                    if m!='diff': cbr = fig.colorbar(sm, ax=axs_list[r*len(method_name):(r+1)*len(method_name)],aspect=27/nrow,location='right')
                else:
                    # will be wonky with m=='diff' choice
                    cbr = fig.colorbar(sm, ax=axs_list,aspect=27/3*nrow)
                    if m=='diff': cbr2 = fig.colorbar(sm2, ax=axs_list,aspect=27/3*nrow)
                if stat == 'kgem': 
                    cbr.ax.set_ylabel(stat_word0)
                else:
                    if do_rel and var!='wallClockTime': 
                        cbr.ax.set_ylabel('relative '+ stat_word0)
                        if m=='diff': cbr2.ax.set_ylabel('relative '+ stat_word2)
                    else:
                        cbr.ax.set_ylabel(stat_word0 + ' [{}]'.format(leg_titl[j]))
                        if m=='diff': cbr2.ax.set_ylabel(stat_word2 + ' [{}]'.format(leg_titl[j]))

                if var=='scalarTotalET' and (stat =='mean' or stat =='mnnz') and m!='diff':
                    # Customizing the tick labels to include negative signs
                    def format_tick(value, tick_number):
                        rounded_value = int(round(value,-2))
                        return f"-{rounded_value}"
                    cbr.ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_tick))
                if m=='diff':
                    # Customizing the tick labels
                    cbr.ax.yaxis.set_major_formatter(ScalarFormatter())

            # lakes
            if m=='ref' and var != 'wallClockTime':
                # only plot wallClockTime for the reference solution
                plt_nm =''
            else:
                if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)

    else: # extra mean/amax variables
        for i,v in enumerate(plot_vars_exVar):
            vmin,vmax = 0, maxes_exVar[i]
            if (v=='airtemp' or v== 'scalarRootZoneTemp' or v=='balanceSoilNrg'): 
                #vmin,vmax = 260, maxes_exVar[i]
                my_cmap2 = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
                my_cmap2.set_bad(color='white') #nan color white
                vmin,vmax =  (273.16-(maxes_exVar[i]-273.16)),maxes_exVar[i],
                norm2 = matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=273.16, vmax=vmax)            
            else:
                norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=0.5)

            r = i//ncol + base_row
            c = i - (r-base_row)*ncol
            m = 'exVar'

            # Plot the data with the full extent of the bas_albers shape
            if (v=='airtemp' or v== 'scalarRootZoneTemp' or v=='balanceSoilNrg'): 
                bas_albers.plot(ax=axs[r,c], column=v+m, edgecolor='none', legend=False, cmap=my_cmap2, norm=norm2,zorder=0)
            else:
                bas_albers.plot(ax=axs[r,c], column=v+m, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)
            stat_word0 = stat_word
            print(f"{'all HRU mean for '}{v+m:<35}{np.nanmean(bas_albers[v+m].values):<10.5f}{' max: '}{np.nanmax(bas_albers[v+m].values):<10.5f}")
            axs[r,c].set_title(plt_name[i])
            axs[r,c].axis('off')
            axs[r,c].set_xlim(xmin, xmax)
            axs[r,c].set_ylim(ymin, ymax)
   
            if (v=='airtemp' or v== 'scalarRootZoneTemp' or v=='balanceSoilNrg'):
                sm = matplotlib.cm.ScalarMappable(cmap=my_cmap2, norm=norm2)
            else:
                sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            sm.set_array([])
            if i==len(plot_vars_exVar)-1: 
                pad = 0.05
            elif i==len(plot_vars_exVar)-2: 
                pad = -0.05
            else: 
                pad = -0.5
            if one_plot:
                cbr = fig.colorbar(sm,ax=axs_list[r*ncol:r*ncol+c+1],aspect=27/nrow, pad=pad)
            else:
                cbr = fig.colorbar(sm,ax=axs_list[r*ncol:r*ncol+c+1],aspect=27/nrow, pad=pad)
            cbr.ax.set_ylabel(stat_word0 + ' [{}]'.format(leg_titl_exVar[i]))

            # lakes
            if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)



# Specify plotting options
if one_plot:
    use_vars = [1]
    use_meth = [0,2,4]
    use_vars_exVar = [3,0,1]
else:
    use_vars = [0,1,2,3,4,5]
    use_vars = [1,5]
    use_meth = [0,1,2,3,4,6]
    use_vars_exVar = [3,0,2,1]
if more_mean: 
    use_vars = ['exVar'] + use_vars # 'exVar' is the extra variables in a balance file, all same method
    if len(use_meth) < len(use_vars_exVar): use_vars_exVar = use_vars_exVar[:len(use_meth)] # chop if longer
    plot_vars_exVar = [plot_vars_exVar[i] for i in use_vars_exVar]
    leg_titl_exVar = [leg_titl_exVar[i] for i in use_vars_exVar]
    maxes_exVar = [maxes_exVar[i] for i in use_vars_exVar]
    #plt_name_exVar = [f"({chr(97+n)}) {plt_titl_exVar[i]+ ' ' + plt_name0_exVar}" for n,i in enumerate(use_vars_exVar)]
    plt_name_exVar = [f"({chr(97+n)}) {plt_titl_exVar[i]}" for n,i in enumerate(use_vars_exVar)]

plot_vars = [var + '_calc' if c == 1 else var for var, c in zip(plot_vars, calc)]
plot_vars = [plot_vars[i] if i != 'exVar' else 'exVar' for i in use_vars]
plt_titl = [plt_titl[i] if i != 'exVar' else 'exVar' for i in use_vars]
leg_titl = [leg_titl[i] if i != 'exVar' else 'exVar' for i in use_vars]
maxes = [maxes[i] if i != 'exVar' else 'exVar' for i in use_vars]
method_name = [method_name[i] for i in use_meth]

if one_plot:
    ncol = len(use_meth)
    nrow = len(use_vars)

    # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 33})
        if more_mean:
            fig,axs = plt.subplots(nrow,ncol,figsize=(16.9*ncol,13*nrow),constrained_layout=True)
        else:
            fig,axs = plt.subplots(nrow,ncol,figsize=(15*ncol,13*nrow),constrained_layout=True)

    else:
        plt.rcParams.update({'font.size': 120})
        if more_mean:
            fig,axs = plt.subplots(nrow,ncol,figsize=(80*ncol,58*nrow),constrained_layout=True)
        else:
            fig,axs = plt.subplots(nrow,ncol,figsize=(67*ncol,58*nrow),constrained_layout=True)

    axs_list = axs.ravel().tolist()
    fig.suptitle('hourly statistics', fontsize=40,y=1.05)
    plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

else:
    #size hardwired to 2x3 for now
    ncol = 3
    nrow = 2
    if len(method_name)>6:
        print('Too many methods for 3x2 plot')
        sys.exit()

    base_row = 0
    plt_name_orig = [f"({chr(97+n)}) {plt_name0[i]}" for n,i in enumerate(use_meth)]

for i,(var,the_max) in enumerate(zip(plot_vars,maxes)):
    
    if one_plot:
        # Reset the names
        base_row = i
        if (len(use_vars)>1): plt_name = [f"({chr(97+n+i*len(use_meth))}) {plt_titl[i] + ' ' + plt_name0[j]}" for n,j in enumerate(use_meth)]
        if (len(use_vars)==1): plt_name = [f"({chr(97+n+i*len(use_meth))}) {plt_name0[j]}" for n,j in enumerate(use_meth)]
        if(var=='exVar'): plt_name = plt_name_exVar
    else:
        plt_name = plt_name_orig
        if(var=='exVar'): plt_name = plt_name_exVar
        # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
        if 'compressed' in fig_fil:
            plt.rcParams.update({'font.size': 33})
            fig,axs = plt.subplots(nrow,ncol,figsize=(15*ncol,13*nrow),constrained_layout=True)
        else:
            plt.rcParams.update({'font.size': 120})
            fig,axs = plt.subplots(nrow,ncol,figsize=(67*ncol,58*nrow),constrained_layout=True)

        # Remove the extra subplots
        if len(method_name) < nrow*ncol:
            for j in range(len(method_name),nrow*ncol):
                r = j//ncol
                c = j-r*ncol
                fig.delaxes(axs[r, c])

        axs_list = axs.ravel().tolist()
        fig.suptitle('{} hourly statistics'.format(plt_titl[i]), fontsize=40,y=1.05)
        plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

    run_loop(i,var,the_max)

    if not one_plot:
        # Save the figure
        fig_fil1 = (var+fig_fil).format(stat)
        plt.savefig(viz_dir/fig_fil1, bbox_inches='tight', transparent=True)

if one_plot:
    # Save the figure
    fig_fil1 = ('all'+fig_fil).format(stat)
    plt.savefig(viz_dir/fig_fil1, bbox_inches='tight', transparent=True)
