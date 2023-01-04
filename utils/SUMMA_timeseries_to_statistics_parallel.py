# SUMMA can produce split-domain output files if the model is used with the -g <start> <num> command line option.
'''Loads timeseries of simulated variables and computes a variety of statistics.'''

# This script analyzes the resulting files and summarizes the timeseries of a (set of) variable(s) into a statistical value.
# Currently implemented are finding the maximum and mean value across the entire time series.
# Outputs are stored in a single file that covers the full spatial extent of the domain.
# Note that this requires the Python package multiprocessing, which is not included in the provided environment.yml and requirements.txt files.
# written originally by W. Knoben, modified by A. Van Beusekom (2023)

import os
import glob
import xarray as xr
from pathlib import Path
import multiprocessing as mp

# Settings
method_name = 'sundials_1en6'
bench_name = 'cat_sundials_1en8'
src_dir = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-' + method_name
ben_dir = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-' + bench_name
src_pat = 'run1_G*_timestep.nc'
des_dir = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/statistics'
des_fil = method_name + '_hrly_diff_stats_{}_{}_{}.nc'
stat = 'max'
settings= {'averageRoutedRunoff': stat; 'wallClockTime': stat; 'scalarTotalET': stat; 'scalarSWE': stat; 'scalarCanopyWat': stat; 'scalarTotalSoilWat': stat}

viz_fil = method_name + '_hrly_diff_stats_{}_{}.nc'
viz_fil = viz_fil.format(','.join(settings.keys()),','.join(settings.values()))

# Make sure we're dealing with the right kind of inputs
src_dir = Path(src_dir)
ben_dir = Path(ben_dir)
des_dir = Path(des_dir)

# Ensure the output path exists
des_dir.mkdir(parents=True, exist_ok=True)

# Get the names of all inputs, assuming folders have same splits of domains and same file names
src_files = glob.glob(str( src_dir / src_pat ))
src_files.sort()

# -- functions
def run_loop(file,bench):

    # extract the subset IDs
    subset = file.split('/')[-1].split('_')[1]

    # open file
    dat = xr.open_dataset(file)
    ben = xr.open_dataset(bench)

    diff = (np.fabs(dat - ben))
    # get rid of gru dimension, assuming they are same as the often are (everything now as hruId)
    diff = diff.drop_vars(['hruId','gruId']) #'wallClockTime'
    m = diff.drop_dims('hru')
    m = m.rename({'gru': 'hru'})
    diff = diff.drop_dims('gru')
    diff = xr.merge([diff,m])

    # compute the requested statistics
    for var,stat in settings.items():

        # Select the case
        if stat == 'mean':
            new = diff[var].mean(dim='time')
        elif stat == 'max':
            new = diff[var].max(dim='time')

        new.to_netcdf(des_dir / des_fil.format(stat,var,subset))
    return

def merge_subsets_into_one(src,pattern,des,name):

    '''Merges all files in {src} that match {pattern} into one file stored in /{des}/{name.nc}'''

    # Find all files
    src_files = glob.glob(str( src / pattern ))

    # Merge into one
    out = xr.merge([xr.open_dataset(file) for file in src_files])

    # save to file
    out.to_netcdf(des / name)

    return #nothing
# -- end functions

# -- start parallel processing
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
if __name__ == "__main__":
    pool = mp.Pool(processes=ncpus)
    pool.map(run_loop,src_files)
    pool.close()
# -- end parallel processing

# merge the individual files into one for further vizualization
merge_subsets_into_one(des_dir,des_fil.replace('{}','*'),des_dir,viz_fil)

# remove the individual files for cleanliness
for file in glob.glob(des_fil.replace('{}','*')):
    os.remove(file)