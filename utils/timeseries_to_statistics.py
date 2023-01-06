# SUMMA can produce split-domain output files if the model is used with the -g <start> <num> command line option.
'''Loads timeseries of simulated variables and computes a variety of statistics.'''

# This script analyzes the resulting files and summarizes the timeseries of a (set of) variable(s) into a statistical value.
# Currently implemented are finding the maximum and mean value across the entire time series.
# Outputs are stored in a single file that covers the full spatial extent of the domain.
# written originally by W. Knoben, modified by A. Van Beusekom (2023)
#   Best to comment out parallel processing lines and run that way on Graham or for full dataset

# Run:
# module load gcc/9.3.0
# module load scipy-stack
# module load gdal/3.2.3
# virtualenv --no-download $SLURM_TMPDIR/env
# source $SLURM_TMPDIR/env/bin/activate
# pip install --no-index --upgrade pip
# pip install --no-index -r requirements.txt
# python timeseries_to_statistics.py


import os
import glob
import xarray as xr
from pathlib import Path
import numpy as np

# Settings
method_name = 'sundials_1en6'
bench_name  = 'sundials_1en8_cat'
top_fold    = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/'

testing = False
if testing: 
    top_fold = '/Users/amedin/Research/USask/test_py/'
else:
    import multiprocessing as mp

src_dir =  top_fold + 'summa-' + method_name
ben_dir =  top_fold + 'summa-' + bench_name
src_pat = 'run1_G*_timestep.nc'
des_dir =  top_fold + 'statistics'
des_fil = method_name + '_hrly_diff_stats_{}_{}_{}.nc'
stat = 'max'
settings= {'averageRoutedRunoff': stat, 'wallClockTime': stat, 'scalarTotalET': stat, 'scalarSWE': stat, 'scalarCanopyWat': stat, 'scalarTotalSoilWat': stat}

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
ben_files = glob.glob(str( ben_dir / src_pat ))
ben_files.sort()

assert len(ben_files) == len(src_files), \
    'Found {} files but need {}!'.format(len(src_files), len(ben_files))


# -- functions
def run_loop(file,bench):
    
    # extract the subset IDs
    subset = file.split('/')[-1].split('_')[1]

    # open file
    dat,ben = xr.open_dataset(file), xr.open_dataset(bench)

    diff = (np.fabs(dat - ben))
    # get rid of gru dimension, assuming they are same as the often are (everything now as hruId)
    diff = diff.drop_vars(['hruId','gruId'])
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
        
    print("wrote output: %s" % (top_fold + 'statistics/' +subset))
        
    return #nothing

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


if testing: 
    # -- no parallel processing
    for (file, bench) in zip(src_files,ben_files):   
        run_loop(file, bench)
else:
    # -- start parallel processing
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    if __name__ == "__main__":
        pool = mp.Pool(processes=ncpus)
        results = [pool.apply_async(run_loop, args=(file, bench)) for (file, bench) in zip(src_files,ben_files)]
        dojob = [p.get() for p in results]
        pool.close()
    # -- end parallel processing

     
# merge the individual files into one for further vizualization
merge_subsets_into_one(des_dir,des_fil.replace('{}','*'),des_dir,viz_fil)

# remove the individual files for cleanliness
for file in glob.glob(str(des_dir / des_fil.replace('{}','*'))):
    os.remove(file)
    
    