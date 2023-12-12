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
# python plot_per_GRU.py sundials_1en6 [stat]
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

# The first input argument specifies the run where the files are
method_name = sys.argv[1] # sys.argv values are strings by default so this is fine (sundials_1en6 or be1)
stat = sys.argv[2]

# Simulation statistics file locations
settings= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
viz_dir = Path('/home/avanb/scratch/statistics')
viz_fil = method_name + '_hrly_diff_stats_{}.nc'
viz_fil = viz_fil.format(','.join(settings))
eff_fil = 'eff_' + method_name + '.txt'
nbatch_hrus = 518 # number of HRUs per batch
use_eff = False # use efficiency in wall clock time
do_rel = True # plot relative to the benchmark simulation
if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Specify variables of interest
plot_vars = settings.copy()
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$mm~y^{-1}$','$s$']
leg_titlm= ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~h^{-1}$','$kg~m^{-2}$','$mm~h^{-1}$','$s$']

if stat == 'rmse': 
    maxes = [2,15,250,0.08,200,10e-3] #[2,15,8e-6,0.08,6e-9,10e-3]
    #maxes = [0.25,2,30,0.01,30,2e-3] #[0.25,2,1e-6,0.01,1e-9,2e-3]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,10e-3]
if stat == 'rmnz': 
    maxes = [2,15,250,0.08,200,10e-3]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,10e-3]
if stat == 'maxe': 
    maxes = [15,25,0.8,2,0.3,0.2] #[15,25,25e-5,2,1e-7,0.2]
    if do_rel: maxes = [0.6,0.02,0.6,0.3,0.6,0.2]
if stat == 'kgem': 
    maxes = [0.9,0.9,0.9,0.9,0.9,10e-3]
if stat == 'mean': 
    maxes = [80,1500,1500,3000,10e-3] #[80,1500,5e-5,8,1e-7,10e-3]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,10e-3]
if stat == 'amax': 
    maxes = [240,1800,3.5,25,7.5,0.2] #[240,1800,1e-3,25,2e-6,0.2]
    if do_rel: maxes = [1.1,1.1,1.1,1.1,1.1,0.2]

fig_fil = method_name + '_hrly_diff_stats_{}_{}_compressed.png'
if do_rel: fig_fil = method_name + '_hrly_diff_stats_{}_{}_rel_compressed.png'
fig_fil = fig_fil.format(','.join(settings),stat)

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
summa = xr.open_dataset(viz_dir/viz_fil)
if use_eff:
    # Read the data from the eff.txt file into a DataFrame
    eff = pd.read_csv(viz_dir/eff_fil, sep=',', header=None, names=['CPU Efficiency', 'Array ID', 'Job Wall-clock time', 'Node Number'])
    # Extract only the values after the ':' character in the 'CPU Efficiency', 'Job Wall-clock time', and 'Node Number' columns
    eff['CPU Efficiency'] = eff['CPU Efficiency'].str.split(':').str[1].astype(float)
    eff['Array ID'] = eff['Array ID'].str.split(':').str[1].astype(int)
    eff['Job Wall-clock time'] = eff['Job Wall-clock time'].str.split(':').str[1].astype(float)
    eff['Node Number'] = eff['Node Number'].str.split(':').str[1].astype(int)

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

    if do_rel: s_rel = summa[plot_var].sel(stat=statr)
    s = summa[plot_var].sel(stat=stat0)
    if do_rel and plot_var != 'wallClockTime': s = s/s_rel

    if plot_var == 'wallClockTime' and use_eff:
        batch = np.floor(np.arange(len(s.indexes['hru'])) /nbatch_hrus)
        #basin_num = np.arange(len(s.indexes['hru'])) % nbatch_hrus #not currently using
        # Create a dictionary to store the values for each batch
        efficiency = {}
        # Iterate over the rows in the data DataFrame
        for index, row in eff.iterrows():
            # Extract the values from the row
            batch0 = int(row['Array ID'])
            eff0 = row['CPU Efficiency']
            # Store the value for the current batch in the dictionary
            efficiency[batch0] = eff0
        # Select the values for the current batch using boolean indexing
        eff_batch = np.array([efficiency[b] for b in batch])
        #node_batch = np.array([node[b] for b in batch]) #not currently using
        # Multiply the s values by efficiency
        s = s*eff_batch

    # Make absolute value norm, not all positive
    s = np.fabs(s) 
    
    # Replace inf values with NaN in the s DataArray
    s = s.where(~np.isinf(s), np.nan)

    if plot_var == 'scalarTotalET' and not do_rel:
        if stat =='rmse' or stat =='rmnz' : s = s*31557600 # make annual total
        if stat =='maxe': s = s*3600 # make hourly max
    if plot_var == 'averageRoutedRunoff' and not do_rel:
        if stat =='rmse' or stat =='rmnz' : s = s*31557600*1000 # make annual total
        if stat =='maxe': s = s*3600*1000 # make hourly max    
    bas_albers[plot_var] = s.sel(hru=hru_ids_shp.values)

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

# Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(3,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(3,2,figsize=(140,133))
fig.suptitle('{} Hourly Statistics'.format(method_name), fontsize=40,y=0.95)

plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

# colorbar axes
f_x_mat = [0.443,0.94,0.443,0.94,0.443,0.94]
f_y_mat = [0.71,0.71,0.38,0.38,0.047,0.047]

plt.tight_layout()

def run_loop(i,var,the_max,f_x,f_y):
    my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
    my_cmap.set_bad(color='white') #nan color white
    vmin,vmax = 0, the_max
    if stat =='mean' and var=='scalarTotalSoilWat' and not do_rel: vmin,vmax = 700, the_max
    if stat =='amax' and var=='scalarTotalSoilWat' and not do_rel: vmin,vmax = 1000, the_max
    if (stat == 'mean' or stat == 'mnnz' or stat == 'amax') and var!='wallClockTime' and do_rel: vmin,vmax = 0.9, the_max
 
    norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=0.5)
    if stat =='kgem' and var!='wallClockTime':
        my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno')) # copy the default cmap
        my_cmap.set_bad(color='white') #nan color white
        vmin,vmax = the_max, 1.0
        norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=1.5)
    r = i//2
    c = i-r*2

    # Plot the data with the full extent of the bas_albers shape
    bas_albers.plot(ax=axs[r,c], column=var, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)

    if stat0 == 'rmse': stat_word = 'RMSE'
    if stat0 == 'rmnz': stat_word = 'RMSE' # no 0s'
    if stat0 == 'maxe': stat_word = 'max abs error'
    if stat0 == 'kgem': stat_word = 'KGE"'
    if stat0 == 'mean': stat_word = 'abs mean'
    if stat0 == 'mnnz': stat_word = 'abs mean' # no 0s'
    if stat0 == 'amax': stat_word = 'abs max'

    if statr == 'mean_ben': statr_word = 'mean'
    if statr == 'mnnz_ben': statr_word = 'mean' # no 0s'
    if statr == 'amax_ben': statr_word = 'max'

    axs[r,c].set_title(plt_titl[i])
    axs[r,c].axis('off')
    axs[r,c].set_xlim(xmin, xmax)
    axs[r,c].set_ylim(ymin, ymax)

    # Custom colorbar
    cax = fig.add_axes([f_x,f_y,0.02,0.25])
    sm = matplotlib.cm.ScalarMappable(cmap=my_cmap, norm=norm)
    sm._A = []
    cbr = fig.colorbar(sm, cax=cax) #, extend='max') #if max extend can't get title right
    if stat == 'rmse' or stat == 'rmnz' or stat == 'mean': cbr.ax.set_ylabel(stat_word + ' [{}]'.format(leg_titl[i]), labelpad=40, rotation=270)
    if stat == 'maxe' or stat == 'amax': cbr.ax.set_ylabel(stat_word + ' [{}]'.format(leg_titlm[i]), labelpad=40, rotation=270)
    if stat == 'kgem': cbr.ax.set_ylabel(stat_word, labelpad=40, rotation=270)
    #if do_rel and var!='wallClockTime': cbr.ax.set_ylabel(stat_word + ' rel to bench ' + statr_word, labelpad=40, rotation=270)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    #cbr.ax.yaxis.set_offset_position('right')

    # lakes
    if plot_lakes: large_lakes_albers.plot(ax=axs[r,c], color=lake_col, zorder=1)

for i,(var,the_max,f_x,f_y) in enumerate(zip(plot_vars,maxes,f_x_mat,f_y_mat)):
    run_loop(i,var,the_max,f_x,f_y)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)
