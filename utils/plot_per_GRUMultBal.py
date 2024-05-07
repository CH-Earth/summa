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
# where stat is rmse or maxe or kgem


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

# plot all runs, pick statistic
stat = sys.argv[1]
method_name=['be1','be1cm','be1en','sundials_1en6cm']  #maybe make this an argument
plt_name=['(a) SUMMA-BE1 Common Form of Heat Eq.','(b) SUMMA-BE1 Temperature Form of Heat Eq.','(c) SUMMA-BE1 Mixed Form of Heat Eq.','(d) SUMMA-SUNDIALS Temperature Form of Heat Eq.']

# Simulation statistics file locations
settings= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']

viz_dir = Path('/home/avanb/scratch/statistics')
viz_fil = method_name.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(['balance']))
do_rel = False # use scaled values
nbatch_hrus = 518 # number of HRUs per batch

# Specify variables of interest
plt_titl = ['Canopy Air Space Energy Balance','Vegetation Energy Balance','Snow Energy Balance','Soil Energy Balance','Vegetation Mass Balance','Snow Mass Balance','Soil Mass Balance','Aquifer Mass Balance', 'Wall Clock Time']
leg_titl = ['$W~m^{-3}$'] * 4 +['$num$'] + ['$kg~m^{-2}~s^{-1}$'] * 4

fig_fil= '_hrly_balance_{}_compressed.png'
if do_rel: 
    fig_fil = '_hrly_scaledBalance_{}_rel_compressed.png'
    for i in range(8):
        settings[i] = 'scaledB' + settings[i][1:]
        plt_titl[i] = 'Scaled ' + plt_titl[i]
    leg_titl = ['$s^{-1}$'] * 8 + ['$s$']

plot_vars = settings.copy()

if stat == 'mean': 
    maxes = [1e-4,1e0,1e0,1e0]+[1e-12,1e-11,1e-10,1e-13] + [3e-3]
    if do_rel: maxes = [1e-6,1e-4,1e-6,1e-7]+[1e-10,1e-11,1e-13,1e-11] + [3e-3]
if stat == 'amax': 
    maxes = [1e-3,1e3,1e3,1e2]+[1e-11,1e-6,1e-7,1e-8] + [1e0]
    if do_rel: maxes = [1e-2,1e0,1e-4,1e-2]+[1e-7,1e-8,1e-10,1e-6] + [1e0]

# Get the albers shapes
main = Path('/home/avanb/projects/rpp-kshook/wknoben/CWARHM_data/domain_NorthAmerica/shapefiles/albers_projection')

# Plot lakes?
plot_lakes = True
# lakes shapefile WHERE IS THIS
#lake_path = Path('C:/Globus endpoint/HydroLAKES/HydroLAKES_polys_v10_shp')
#lake_name = 'HydroLAKES_polys_v10_subset_NA.shp'

## Control file handling
# Store the name of the 'active' file in a variable
controlFile = 'plot_control_NorthAmerica.txt'

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

## Catchment shapefile location and variable names
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
# Plot rivers?
plot_rivers = False
# River network path & name
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
# Set the target CRS
acc = 'ESRI:102008'

# catchment shapefile, first 2 lines throw error so cutting them
#bas = gpd.read_file(hm_catchment_path/hm_catchment_name)
#bas_albers = bas.to_crs(acc)
bas_albers = gpd.read_file(main/'basin.shp')
xmin, ymin, xmax, ymax = bas_albers.total_bounds

# river network shapefile, first 2 lines throw error so cutting them
if plot_rivers:
	#riv = gpd.read_file(river_network_path/river_network_name)
	#riv_albers = riv.to_crs(acc)
	riv_albers = gpd.read_file(main/'river.shp')

# lakes shapefile, first 2 lines throw error so cutting them
if plot_lakes:
	#lakes = gpd.read_file(lake_path/lake_name)
	#lak_albers = lakes.to_crs(acc)
	lak_albers = gpd.read_file(main/'lakes.shp')

## Pre-processing, map SUMMA sims to catchment shapes
# Get the aggregated statistics of SUMMA simulations
summa = {}
for i, m in enumerate(method_name):
    # Get the aggregated statistics of SUMMA simulations
    summa[m] = xr.open_dataset(viz_dir/viz_fil[i])

# Match the accummulated values to the correct HRU IDs in the shapefile
hru_ids_shp = bas_albers[hm_hruid].astype(int) # hru order in shapefile
for plot_var in plot_vars:
    stat0 = stat

    for m in method_name:
        s = summa[m][plot_var].sel(stat=stat0)

        # Make absolute value norm, not all positive
        s = np.fabs(s) 

        # Replace inf and 9999 values with NaN in the s DataArray
        s = s.where(~np.isinf(s), np.nan).where(lambda x: x != 9999, np.nan)
        
        # Create a new column in the shapefile for each method, and fill it with the statistics
        bas_albers[plot_var+m] = np.nan
        hru_ind = [i for i, hru_id in enumerate(hru_ids_shp.values) if hru_id in s.hru.values] # if some missing
        bas_albers.loc[hru_ind, plot_var+m] = s.sel(hru=hru_ids_shp.values[hru_ind]).values 
        #bas_albers[plot_var+m]= s.sel(hru=hru_ids_shp.values)

# Select lakes of a certain size for plotting
if plot_lakes:
    minSize = 1000 # km2
    in_domain = (lak_albers['Country'] == 'Canada') | \
                (lak_albers['Country'] == 'United States of America') | \
                (lak_albers['Country'] == 'Mexico')
    out_domain = (lak_albers['Pour_long'] > -80) & (lak_albers['Pour_lat'] > 65) # Exclude Baffin Island
    large_lakes_albers = lak_albers.loc[(lak_albers['Lake_area'] > minSize) & in_domain & (~out_domain) ]
    lake_col = (8/255,81/255,156/255)

##Figure

def run_loop(j,var,the_max):
    stat0 = stat

    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    my_cmap.set_bad(color='white') #nan color white
    vmin,vmax = the_max*1e-4, the_max
    if any(substring in var for substring in ['VegNrg', 'SnowNrg', 'SoilNrg']):
        vmin, vmax = the_max * 1e-7, the_max
    if var in ['wallClockTime',]: vmin,vmax = the_max*1e-1, the_max
 
    norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'amax': stat_word = 'max'

    # colorbar axes
    f_x_mat = [0.46,0.96,0.46,0.96]
    f_y_mat = [0.55,0.55,0.07,0.07]

    for i,(m,f_x,f_y) in enumerate(zip(method_name,f_x_mat,f_y_mat)):
        r = i//2
        c = i-r*2

        # Plot the data with the full extent of the bas_albers shape
        bas_albers.plot(ax=axs[r,c], column=var+m, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)

        axs[r,c].set_title(plt_name[i])
        axs[r,c].axis('off')
        axs[r,c].set_xlim(xmin, xmax)
        axs[r,c].set_ylim(ymin, ymax)

        # Custom colorbar
        cax = fig.add_axes([f_x,f_y,0.02,0.375])
        sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
        sm.set_array([])
        cbr = fig.colorbar(sm, cax=cax) #, extend='max') #if max extend can't get title right
        cbr.ax.set_ylabel(stat_word + ' [{}]'.format(leg_titl[j]), labelpad=40, rotation=270)
        if do_rel and var!='wallClockTime': cbr.ax.set_ylabel('scaled '+ stat_word, labelpad=40, rotation=270)

        #cbr.ax.yaxis.set_offset_position('right')

        # lakes
        if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)

for i,(var,the_max) in enumerate(zip(plot_vars,maxes)):
 
    # Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
    if 'compressed' in fig_fil:
        plt.rcParams.update({'font.size': 25})
    else:
        plt.rcParams.update({'font.size': 100})

    if 'compressed' in fig_fil:
        fig,axs = plt.subplots(2,2,figsize=(35,28))
    else:
        fig,axs = plt.subplots(2,2,figsize=(140,133))

    # Remove the fourth subplot
    fig.delaxes(axs[1, 1])

    fig.suptitle('{} Hourly Statistics'.format(plt_titl[i]), fontsize=40,y=1.05)

    plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

    plt.tight_layout()
 
    run_loop(i,var,the_max)

    fig_fil1 = (var+fig_fil).format(stat)
    # Save
    plt.savefig(viz_dir/fig_fil1, bbox_inches='tight', transparent=True)
