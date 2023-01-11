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
# python plot_per_GRU.py sundials_1en6 rmse


# modules
import sys
import os
import pyproj 
import fiona
import matplotlib
import numpy as np
import xarray as xr
import geopandas as gpd
from pathlib import Path
import matplotlib.pyplot as plt
import copy
# for parallellization
import multiprocessing as mp
from functools import partial
import io
from copy import deepcopy
from PIL import Image

# The first input argument specifies the run where the files are
method_name = sys.argv[1] # sys.argv values are strings by default so this is fine (sundials_1en6 or be1)
stat = sys.argv[2] # max or rmse

# Simulation statistics file locations
settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}
#viz_dir = Path('/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics')
viz_dir = Path('/home/avanb/scratch/statistics')
viz_fil = method_name + '_hrly_diff_stats_{}_{}.nc'
viz_fil = viz_fil.format(','.join(settings.keys()),','.join(settings.values()))

# Specify variables of interest
plot_vars = ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','averageRoutedRunoff','wallClockTime']
plt_titl = ['(a) Snow Water Equivalent','(b) Total soil water content','(c) Total evapotranspiration', '(d) Total water on the vegetation canopy','(e) Average routed runoff','(f) Wall clock time']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','$kg~m^{-2}~s^{-1}$','$kg~m^{-2}$','$m~s^{-1}$','$s$']
if stat=='rmse': maxes = [60,30,5e-6,0.2,1e-8,200]
if stat=='max' : maxes = [100,50,5e-4,4,35e-8,50000]

fig_fil = method_name + '_hrly_diff_stats_{}_{}_compressed2.png'
fig_fil = fig_fil.format(','.join(settings.keys()),','.join(settings.values()))

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

# Match the accummulated values to the correct HRU IDs in the shapefile
hru_ids_shp = bas_albers[hm_hruid].astype(int) # hru order in shapefile
for plot_var in plot_vars:
    bas_albers[plot_var] = summa[plot_var].sel(hru=hru_ids_shp.values)

# Select lakes of a certain size for plotting
if plot_lakes:
    minSize = 1000 # km2
    in_domain = (lak_albers['Country'] == 'Canada') | \
                (lak_albers['Country'] == 'United States of America') | \
                (lak_albers['Country'] == 'Mexico')
    out_domain = (lak_albers['Pour_long'] > -80) & (lak_albers['Pour_lat'] > 65) # Exclude Baffin Island
    large_lakes_albers = lak_albers.loc[(lak_albers['Lake_area'] > minSize) & in_domain & (~out_domain) ]
# Set the lake color
if plot_lakes:
    lake_col = (8/255,81/255,156/255)

##Figure

# Set the font size: we need this to be huge so we can also make our plotting area huge, to avoid a gnarly plotting bug
if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 25})
else:
    plt.rcParams.update({'font.size': 100})

# Flip the evaporation values so that they become positive, not if plotting diffs
#bas_albers['plot_ET'] = bas_albers['scalarTotalET'] * -1
#bas_albers['plot_ET'] = bas_albers['plot_ET'].where(bas_albers['scalarTotalET'] != -9999, np.nan)

plt.rcParams['patch.antialiased'] = False # Prevents an issue with plotting distortion along the 0 degree latitude and longitude lines

# colorbar axes
f_x_mat = [0.473,0.97,0.473,0.97,0.473,0.97]
f_y_mat = [0.69,0.69,0.36,0.36,0.03,0.03]

plt.tight_layout()

my_cmap = copy.copy(matplotlib.cm.get_cmap('inferno_r')) # copy the default cmap
#my_cmap.set_bad(my_cmap.colors[0])

if stat == 'rmse':
    stat_word = ' Hourly RMSE'
if stat == 'max':
    stat_word = ' Hourly max abs error'
    
def make_plot(var,the_max,f_x,f_y,plt_t,leg_t, fig, axes):
    #vmin = (bas_albers[var].where(bas_albers[var]>0)).min()
    #a = bas_albers[var].where(bas_albers[var]>0)
    #bas_albers[var] = a.fillna(vmin)
    #vmin,vmax = bas_albers[var].min(), bas_albers[var].max()
    vmin,vmax = 0, the_max
    norm=matplotlib.colors.PowerNorm(vmin=vmin,vmax=vmax,gamma=0.5)
    
    # Data
    sm = bas_albers.plot(ax=axes,column=var, edgecolor='none', legend=False, cmap=my_cmap, norm=norm,zorder=0)
    
    # Custom colorbar
    cax = fig.add_axes([f_x,f_y,0.02,0.3])
    cbr = fig.colorbar(sm, cax=cax, extend='max')
    cbr.ax.set_title('[{}]'.format(leg_t))
    
    # wall Clock doesn't do difference
    if var == 'wallClockTime':
        if stat == 'rmse': stat_word = ' Hourly mean'
        if stat == 'rmse': stat_word = ' Hourly max'

    axes.set_title(plt_t+ stat_word)
    axes.axis('off')
    
    # lakes
    if plot_lakes: large_lakes_albers.plot(color=lake_col, zorder=1)


#From https://github.com/paulgavrikov/parallel-matplotlib-grid
def parallel_plot(plot_fn, data1,data2,data3,data4,data5,data6):
    if 'compressed' in fig_fil:
        fig,axes = plt.subplots(3,2,figsize=(35,33))
    else:
        fig,axes = plt.subplots(3,2,figsize=(140,133))

    
    with mp.Pool() as pool:
        for ax, rastered in zip(axes.ravel(), pool.map(partial(_parallel_plot_worker, plot_fn=plot_fn), data1,data2,data3,data4,data5,data6)):
            ax.imshow(rastered)
            
            # The following code hides axes
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
    
    plt.subplots_adjust(hspace=0, wspace=0)
    
def _parallel_plot_worker( data1,data2,data3,data4,data5,data6, plot_fn):
    fig = plt.figure()
    matplotlib.font_manager._get_font.cache_clear()  # necessary to reduce text corruption artifacts
    axes = plt.axes()
    
    plot_fn( data1,data2,data3,data4,data5,data6, fig, axes)
    pil_img = rasterize(fig)
    plt.close()
    
    return pil_img


def rasterize(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    pil_img = deepcopy(Image.open(buf))
    buf.close()
    
    return pil_img

# Run
parallel_plot(make_plot,plot_vars,maxes,f_x_mat,f_y_mat,plt_titl,leg_titl)


# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)
