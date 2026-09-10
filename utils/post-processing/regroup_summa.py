# regroup summa runs output into different group size
# written by A. Van Beusekom (2025)

# Run:
# python regroup_summa.py sundials_1en8
# or python regroup_summa.py sundials_1en8 $SLURM_ARRAY_TASK_ID

import os
from glob import glob
import netCDF4 as nc
import numpy as np
import sys

# number of groups to divide the output into
num_groups = 65
gru_num = 517315  # Number of GRUs 

missing = False  # if appending nan hrus to batch because failed
missgru = 72055933  # batch 205 summa-be32 value
misshru = missgru  # could be different

run_batch = True # run by batch
if run_batch:
    top_fold = '/project/k/kshook/avanb/enthalpy_paper/runs/'
    job_index = int(sys.argv[2])  # Job array index
else: # run with python parallel processing
    import multiprocessing as mp
    top_fold = '/project/k/kshook/avanb/enthalpy_paper/runs/'

method_name = sys.argv[1]  # sys.argv values are strings by default so this is fine (sundials_1en8 or be64)
ncdir = top_fold + 'summa-' + method_name + '_nocat'
file_pattern = 'run1_G*_timestep.nc'
ctdir = top_fold + 'summa-' + method_name

# Get list of split summa output files (hardwired pattern)
outfilelist0 = glob((ncdir + '/' + file_pattern))
outfilelist0.sort()

# Variables to exclude
exclude_vars = [
    'scalarCanopyWat', 'scalarSWE', 'balanceCasNrg', 'balanceVegNrg', 
    'balanceSnowNrg', 'balanceVegMass', 'balanceSnowMass', 'balanceSoilMass', 
    'balanceAqMass', 'scalarTotalET'
]

# -- functions
def concatenate_files_in_range(outfilelist0, ctdir, start_gru, end_gru):
    out_name = f'run1__G{start_gru:06d}-{end_gru:06d}_timestep.nc'
    filtered_files = [file for file in outfilelist0 if int(file.split('/')[-1].split('_')[1][1:7]) <= end_gru and int(file.split('/')[-1].split('_')[1][8:14]) >= start_gru]
    gru_num = end_gru - start_gru + 1
    hru_num = end_gru - start_gru + 1

    # Write output
    with nc.Dataset(filtered_files[0]) as src:
        with nc.Dataset(ctdir + '/' + out_name, "w") as dst:
            # copy dimensions
            for name, dimension in src.dimensions.items():
                if name == 'gru':
                    dst.createDimension(name, gru_num)
                elif name == 'hru':
                    dst.createDimension(name, hru_num)
                else:
                    dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            # copy variable attributes all at once via dictionary
            gru_vars = [] # variable name, gru axis in variable dimension for concatenation.
            hru_vars = []
            for name, variable in src.variables.items():
                if name in exclude_vars:
                    continue
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
    
                # Assign different values depending on dimension
                dims = variable.dimensions
                if 'gru' in dims:
                    gru_vars.append([name,dims.index('gru')])
                elif 'hru' in dims:
                    hru_vars.append([name,dims.index('hru')])
                else:
                    dst[name][:]=src[name][:]

            # read values for gru and hru dimensioned variables
            Dict = {}
            gru_vars_num = len(gru_vars)
            hru_vars_num = len(hru_vars)
            for i, file in enumerate(filtered_files):
                start_file = int(file.split('/')[-1].split('_')[1][1:7])
                end_file = int(file.split('/')[-1].split('_')[1][8:14])
                start = 0
                end = None
                if i==0: start = start_gru-start_file
                if i==len(filtered_files)-1: end = end_gru - end_file
                if end == 0: end = None

                print("combining file %d %s" % (i,file))
                f = nc.Dataset(file)
                for j in range(gru_vars_num):
                    gru_var_name = gru_vars[j][0]
                    dim_index = gru_vars[j][1]
                    data=f[gru_var_name][:]
                    slices = [slice(None)] * data.ndim
                    slices[dim_index] = slice(start, end)
                    data = data[tuple(slices)]
                    if i == 0:
                        Dict[gru_var_name]=data
                    else:
                        Dict[gru_var_name]=np.concatenate((Dict[gru_var_name],data),axis=dim_index)

                for j in range(hru_vars_num):
                    hru_var_name = hru_vars[j][0]
                    dim_index = hru_vars[j][1]
                    data=f[hru_var_name][:]
                    slices = [slice(None)] * data.ndim
                    slices[dim_index] = slice(start, end)
                    data = data[tuple(slices)]
                    if i == 0:
                        Dict[hru_var_name]=data
                    else:
                        Dict[hru_var_name]=np.concatenate((Dict[hru_var_name],data),axis=dim_index)

            # assign values for gru and hru dimensioned variables
            for j in range(gru_vars_num):
                dst.variables[gru_vars[j][0]][:] = Dict[gru_vars[j][0]]
            for j in range(hru_vars_num):
                dst.variables[hru_vars[j][0]][:] = Dict[hru_vars[j][0]]

            #if missing HRUs, this is slow or broken
            if missing:
                new_index = np.append(dst["gru"][:],missgru)
                dst.reindex({"gru": new_index})
                dst.sel(gru=missgru)["gruId"] = missgru

                new_index = np.append(dst["hru"][:],misshru)
                dst.reindex({"hru": new_index})
                dst.sel(gru=misshru)["hruId"] = misshru

        print("wrote output: %s" % (ctdir + '/' + out_name))

    return out_name

# -- end functions

def process_group(g):
    size_g = int(np.ceil(gru_num / num_groups))
    start_gru = g * size_g + 1
    end_gru = min((g + 1) * size_g, gru_num)
    concatenate_files_in_range(outfilelist0, ctdir, start_gru, end_gru)

if run_batch:
    # -- no parallel processing
    process_group(job_index)
else:
    # -- start parallel processing
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    if __name__ == "__main__":
        pool = mp.Pool(processes=ncpus)
        results = [pool.apply_async(process_group, args=(g,)) for g in range(num_groups)]
        dojob = [p.get() for p in results]
        pool.close()