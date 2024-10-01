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
# where stat is mean or amax


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


one_plot = False # true is one plot, false is multiple plots (one per variable)
run_local = False # true is run on local machine (only does testing), false is run on cluster
fixed_Mass_units = False # true is convert mass balance units to kg m-2 s-1, if ran new code with depth in calculation


if run_local: 
    stat = 'mean'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')


method_name=['be1','be1cm','be1en','sundials_1en6cm','sundials_1en6en','sundials_1en8cm']  #maybe make this an argument
plt_name0=['SUMMA-BE1 common thermo. eq.','SUMMA-BE1 temperature thermo. eq.','SUMMA-BE1 mixed thermo. eq.','SUMMA-SUNDIALS temperature thermo. eq.','SUMMA-SUNDIALS enthalpy thermo. eq.','reference solution']

# Simulation statistics file locations
settings= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']

viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(['balance']))
nbatch_hrus = 518 # number of HRUs per batch

# Specify variables in files
plt_titl = ['canopy air space enthalpy balance','vegetation enthalpy balance','snow enthalpy balance','soil enthalpy balance','vegetation mass balance','snow mass balance','soil mass balance','aquifer mass balance', 'wall clock time']
leg_titl = ['$W~m^{-3}$'] * 4 + ['$kg~m^{-2}~s^{-1}$'] * 4 + ['$s$']
if fixed_Mass_units: leg_titl = ['$W~m^{-3}$'] * 4 + ['s^{-1}$'] * 3 + ['m~s^{-1}$'] + ['$s$']

fig_fil= '_hrly_balance_{}_compressed.png'
plot_vars = settings.copy()

if stat == 'mean': 
    maxes = [1e-1,1e1,1e1,1e1]+[1e-7,1e-7,1e-7,1e-9] + [20e-3]
if stat == 'amax':
    maxes = [1e1,1e3,1e3,1e3]+[1e-5,1e-5,1e-5,1e-7] + [2.0]

# Get simulation statistics
summa = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fil[i])



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

    for m in method_name:
        s = summa[m][plot_var].sel(stat=stat0)
        if fixed_Mass_units and 'Mass' in plot_var: s = s/1000 # /density for mass balance

        # Make absolute value norm, not all positive
        s = np.fabs(s) 

        # Replace inf and 9999 values with NaN in the s DataArray
        s = s.where(~np.isinf(s), np.nan).where(lambda x: x != 9999, np.nan)
        
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

    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    my_cmap.set_bad(color='white') #nan color white
    vmin,vmax = the_max*1e-9, the_max
    if any(substring in var for substring in ['VegNrg', 'SnowNrg', 'SoilNrg']):
        vmin, vmax = the_max * 1e-9, the_max
    if var in ['wallClockTime',]: vmin,vmax = the_max*1e-1, the_max
    if fixed_Mass_units and 'Mass' in var: # / density for mass balance
        vmin = vmin/1000
        vmax = vmax/1000
 
    norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'amax': stat_word = 'max'

    for i,m in enumerate(method_name):
        r = i//ncol + base_row
        c = i - (r-base_row)*ncol

        # Plot the data with the full extent of the bas_albers shape
        bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)
        print(f"{'all HRU mean for '}{var+m:<35}{np.nanmean(bas_albers[var+m].values):<10.5f}{' max: '}{np.nanmax(bas_albers[var+m].values):<10.5f}")
        axs[r,c].set_title(plt_name[i])
        axs[r,c].axis('off')
        axs[r,c].set_xlim(xmin, xmax)
        axs[r,c].set_ylim(ymin, ymax)

        # Custom colorbar
        # Custom colorbar
        if i==len(method_name)-1:
            sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
            sm.set_array([])
            if one_plot:
                cbr = fig.colorbar(sm, ax=axs_list[r*len(method_name):(r+1)*len(method_name)],aspect=27/nrow)
            else:
                cbr = fig.colorbar(sm, ax=axs_list,aspect=27/3*nrow)
            cbr.ax.set_ylabel(stat_word + ' [{}]'.format(leg_titl[j]))

        # lakes
        if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)



# Specify plotting options
if one_plot:
    use_vars = [1,2,3]
    use_meth = [0,2,4]
else:
    use_vars = [0,1,2,3,4,5,6,7]
    use_vars = [3]
    use_meth = [0,1,2,3,4,5]
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
    nrow = 3
    if len(method_name)>6:
        print('Too many methods for 3x2 plot')
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
