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


do_rel = False # true is plot relative to the benchmark simulation
one_plot = True # true is one plot, false is multiple plots (one per variable)
run_local = False # true is run on local machine (only does testing), false is run on cluster

if run_local: 
    stat = 'mnnz'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')


# NOTE: method_name 'ref' will plot the reference solution, 'diff' will plot the difference between two simulations
#       method_name 'diff' requires the specification of the two simulations to subtract in the variables from_meth and sub_meth

method_name=['be1','be16','be32','sundials_1en6','ref']
plt_name0=['SUMMA-BE1','SUMMA-BE16','SUMMA-BE32','SUMMA-SUNDIALS','reference solution']
#method_name=['be1','be1cm','be1en','sundials_1en6cm','diff']
#plt_name0=['SUMMA-BE1 Common Form of Heat Eq.','SUMMA-BE1 Temperature Form of Heat Eq.','SUMMA-BE1 Mixed Form of Heat Eq.','SUMMA-SUNDIALS Temperature Form of Heat Eq.','SUMMA-BE1 Mixed Form - Temperature Form']
from_meth = 'be1en' # name of the first simulation in the difference simulation, only used if a method_name is 'diff'
sub_meth = 'be1' # name of the simulation to subtract in the difference simulation, only used if a method_name is 'diff'

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']

viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings))
nbatch_hrus = 518 # number of HRUs per batch
if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Specify variables in files
plot_vars = settings.copy()
plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','average routed runoff','wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$s$']

fig_fil= '_hrly_diff_stats_{}_compressed.png'
if do_rel: fig_fil = '_hrly_diff_stats_{}_rel_compressed.png'

if stat == 'rmse' or stat == 'rmnz': 
    maxes = [2,15,250,0.08,200,10e-3] 
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,10e-3]
if stat == 'maxe': 
    maxes = [15,25,0.8,2,0.3,0.2] #[15,25,25e-5,2,1e-7,0.2]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,0.2]
if stat == 'kgem': 
    maxes = [0.9,0.9,0.9,0.9,0.9,10e-3]
if stat == 'mean' or stat == 'mnnz': 
    maxes = [80,1700,2000,8,5000,10e-3] #[80,1500,5e-5,8,1e-7,10e-3]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,10e-3]
if stat == 'amax': 
    maxes = [240,1800,3.5,25,7.5,0.2] #[240,1800,1e-3,25,2e-6,0.2]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,0.2]

# Get simulation statistics
summa = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    if m!='diff' and m!='ref': summa[m] = xr.open_dataset(viz_dir/viz_fil[i])



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
        'hm_hruid': s.hru.values[range(100)],  # Example HRU IDs
        'geometry': [Point(x, y) for x, y in zip(range(100), range(100))]  # Simple geometries
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
for plot_var in plot_vars:
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
    for m in method_name:
        if m=='diff': 
            s = np.fabs(summa[from_meth][plot_var].sel(stat=stat0)) - np.fabs(summa[sub_meth][plot_var].sel(stat=stat0))
        elif m=='ref':
            s = np.fabs(summa[method_name[0]][plot_var].sel(stat=statr))
        else:
            s = np.fabs(summa[m][plot_var].sel(stat=stat0))
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
        bas_albers[plot_var+m] = np.nan
        hru_ind = [i for i, hru_id in enumerate(hru_ids_shp.values) if hru_id in s.hru.values] # if some missing
        bas_albers.loc[hru_ind, plot_var+m] = s.sel(hru=hru_ids_shp.values[hru_ind]).values 

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

    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    my_cmap.set_bad(color='white') #nan color white
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
    vmin,vmax = -the_max/10, the_max/4
    norm2 = matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
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
    
    for i,m in enumerate(method_name):
        r = i//ncol + base_row
        c = i - (r-base_row)*ncol

        # Plot the data with the full extent of the bas_albers shape
        if m=='diff':
            bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap2, norm=norm2,zorder=0)
            stat_word0 = stat_word+' difference'
        else:
            bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)
            stat_word0 = stat_word

        axs[r,c].set_title(plt_name[i])
        axs[r,c].axis('off')
        axs[r,c].set_xlim(xmin, xmax)
        axs[r,c].set_ylim(ymin, ymax)

        # Custom colorbar
        if i==len(method_name)-1:
            if m=='diff':
                sm = matplotlib.cm.ScalarMappable(cmap=my_cmap2, norm=norm2)
            else:
                sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            sm.set_array([])
            if one_plot:
                cbr = fig.colorbar(sm, ax=axs_list[r*len(method_name):(r+1)*len(method_name)],aspect=27/3)
            else:
                cbr = fig.colorbar(sm, ax=axs_list,aspect=27/3*nrow)
            if stat == 'kgem': 
                cbr.ax.set_ylabel(stat_word0)
            else:
                if do_rel and var!='wallClockTime': 
                    cbr.ax.set_ylabel('relative '+ stat_word0)
                else:
                    cbr.ax.set_ylabel(stat_word0 + ' [{}]'.format(leg_titl[j]))

            if var=='scalarTotalET' and (stat =='mean' or stat =='mnnz'):
                # Customizing the tick labels to include negative signs
                def format_tick(value, tick_number):
                    rounded_value = int(round(value,-2))
                    return f"-{rounded_value}"
                cbr.ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_tick))

        # lakes
        if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)



# Specify plotting options
if one_plot:
    use_vars = [1,2,4]
    use_meth = [0,3,4]
else:
    use_vars = [0,1,2,3,4,5]
    use_meth = [0,1,2,3]
plot_vars = [plot_vars[i] for i in use_vars]
plt_titl = [plt_titl[i] for i in use_vars]
leg_titl = [leg_titl[i] for i in use_vars]
maxes = [maxes[i] for i in use_vars]
method_name = [method_name[i] for i in use_meth]

if one_plot:
    ncol = len(use_meth)
    nrow = len(use_vars)

    # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 33})
        fig,axs = plt.subplots(nrow,ncol,figsize=(15*ncol,13*nrow),constrained_layout=True)
    else:
        plt.rcParams.update({'font.size': 120})
        fig,axs = plt.subplots(nrow,ncol,figsize=(67*ncol,58*nrow),constrained_layout=True)

    axs_list = axs.ravel().tolist()
    fig.suptitle('hourly statistics', fontsize=40,y=1.05)
    plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

else:
    #size hardwired to 2x2 for now
    ncol = 2
    nrow = 2
    if len(method_name)>4:
        print('Too many methods for 2x2 plot')
        sys.exit()

    base_row = 0
    plt_name = [f"({chr(97+n)}) {plt_name0[i]}" for n,i in enumerate(use_meth)]

for i,(var,the_max) in enumerate(zip(plot_vars,maxes)):
    
    if one_plot:
        # Reset the names
        base_row = i
        plt_name = [f"({chr(97+n+i*len(use_meth))}) {plt_titl[i] + ' ' + plt_name0[j]}" for n,j in enumerate(use_meth)]
    else:
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
