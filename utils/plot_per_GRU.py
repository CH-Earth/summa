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
# module load gcc/9.3.0
# module load scipy-stack
# module load gdal/3.2.3
# virtualenv --no-download $SLURM_TMPDIR/env
# source $SLURM_TMPDIR/env/bin/activate
# pip install --no-index --upgrade pip
# pip install --no-index -r requirements.txt
# python plot_per_GRU.py


# modules
import pyproj
import matplotlib
import numpy as np
import xarray as xr
import geopandas as gpd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Simulation statistics file locations
method_name = 'sundials_1en6'
stat = 'max'
settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}
viz_dir = Path('/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics')
viz_fil = method_name + '_hrly_diff_stats_{}_{}.nc'
viz_fil = viz_fil.format(','.join(settings.keys()),','.join(settings.values()))

# Specify variables of interest
plot_vars = ['averageRoutedRunoff','wallClockTime','scalarTotalET','scalarSWE','scalarCanopyWat','scalarTotalSoilWat']
fig_fil = method_name + '_hrly_diff_stats_{}_{}_compressed.png'
fig_fil = fig_fil.format(','.join(settings.keys()),','.join(settings.values()))

# Get the albers shapes
main = Path('/home/avanb/projects/rpp-kshook/wknoben/CWARHM_data/domain_NorthAmerica/shapefiles/albers_projection')

# Plot lakes?
plot_lakes = False
# lakes shapefile WHERE IS THIS
lake_path = Path('C:/Globus endpoint/HydroLAKES/HydroLAKES_polys_v10_shp')
lake_name = 'HydroLAKES_polys_v10_subset_NA.shp'

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

# catchment shapefile
bas = gpd.read_file(hm_catchment_path/hm_catchment_name)
bas_albers = bas.to_crs(acc)
bas_albers = gpd.read_file(main/'basin.shp')

# river network shapefile
if plot_rivers:
	riv = gpd.read_file(river_network_path/river_network_name)
	riv_albers = riv.to_crs(acc)
	riv_albers = gpd.read_file(main/'river.shp')

# lakes shapefile
if plot_lakes:
	lakes = gpd.read_file(lake_path/lake_name)
	lakes.to_crs(acc)
	lak_albers = gpd.read_file(main/'lakes.shp')


# Print the median basin size for curiousity
print('median area = {} m^2'.format(bas['HRU_area'].median() / 10**6))
print('mean area   = {} m^2'.format(bas['HRU_area'].mean() / 10**6))

#median area = 33.06877343600296 m^2
#mean area   = 40.19396140285971 m^2

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

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(3,2,figsize=(35,33))
else:
    fig,axs = plt.subplots(3,2,figsize=(140,133))

# colorbar axes
cax1 = fig.add_axes([0.473,0.60,0.02,0.3])
cax2 = fig.add_axes([0.97 ,0.60,0.02,0.3])
cax3 = fig.add_axes([0.473,0.10,0.02,0.3])
cax4 = fig.add_axes([0.97 ,0.10,0.02,0.3])

plt.tight_layout()

# add maps
var = 'scalarSWE'
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[0,0], column=var, edgecolor='none', legend=True,\
                cmap='Greys_r', cax=cax1, norm=norm, zorder=0)
axs[0,0].set_title('(a) Snow Water Equivalent Absolute '+stat+ ' Diffs')
axs[0,0].axis('off')
cax1.set_ylabel('scalarSWE $[kg~m^{-2}]$',labelpad=-600)

# SM
var = 'scalarTotalSoilWat'
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[0,1], column=var, edgecolor='none', legend=True,\
                cmap='cividis_r', cax=cax2, norm=norm, zorder=0)
axs[0,1].set_title('(b) Total soil water content Absolute '+stat+ ' Diffs')
axs[0,1].axis('off')
cax2.set_ylabel('scalarTotalSoilWat $[kg~m^{-2}]$',labelpad=-600)

# ET
var = 'scalarTotalET'
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[1,0], column=var, edgecolor='none', legend=True,\
                cmap='viridis', cax=cax3, norm=norm, zorder=0)
axs[1,0].set_title('(c) Total evapotranspiration Absolute '+stat+ ' Diffs')
axs[1,0].axis('off')
cax3.set_ylabel('scalarTotalET $[kg~m^{-2}~s^{-1}]$',labelpad=-600)

# CanWat
var = 'scalarCanopyWat'
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[1,0], column=var, edgecolor='none', legend=True,\
                cmap='viridis_r', cax=cax3, norm=norm, zorder=0)
axs[1,1].set_title('(d) Total water on the vegetation canopy Absolute '+stat+ ' Diffs')
axs[1,1].axis('off')
cax3.set_ylabel('scalarCanopyWat $[kg~m^{-2}]$',labelpad=-600)

# Runoff
var = 'averageRoutedRunoff'
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[2,0], column=var, edgecolor='none', legend=True,\
                cmap='Blues', cax=cax3, norm=norm, zorder=0)
axs[2,0].set_title('(e) Routed runoff Absolute '+stat+ ' Diffs')
axs[2,0].axis('off')
cax3.set_ylabel('averageRoutedRunoff $[m~s^{-1}]$',labelpad=-600)

# Clock time
var = 'wallClockTime('
norm = matplotlib.colors.LogNorm(vmin=bas_albers[var].min(), vmax=bas_albers[var].max())
bas_albers.plot(ax=axs[2,1], column=var, edgecolor='none', legend=True,\
                cmap='Greys', cax=cax3, norm=norm, zorder=0)
axs[2,1].set_title('(f) Wall clock time Absolute '+stat+ ' Diffs')
axs[2,1].axis('off')
cax3.set_ylabel('wallClockTime( $[s]$',labelpad=-600)


# lakes
if plot_lakes:
    large_lakes_albers.plot(ax=axs[0,0], color=lake_col, zorder=1)
    large_lakes_albers.plot(ax=axs[0,1], color=lake_col, zorder=1)
    large_lakes_albers.plot(ax=axs[1,0], color=lake_col, zorder=1)
    large_lakes_albers.plot(ax=axs[1,1], color=lake_col, zorder=1)

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=True)
