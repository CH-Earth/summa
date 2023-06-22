# concatenate the outputs of a split domain summa run into fewer groups
# written originally by Manab Saharia, updated by Hongli Liu and Andy Wood, and W. Knoben
# modified by A. Van Beusekom (2023)
#   Best to comment out parallel processing lines and run that way on Graham or for full dataset

# Run:
# python concat_groups_split_summa.py sundials_1en8

import os
from glob import glob
import netCDF4 as nc
import numpy as np

catby_num   = 2 #number of files to cat into one, if had to divide runs from regular batches into sub-batches to finish in 7 days
top_fold    = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/'

missing = False # if appending nan hrus to batch because failed
missgru = 72055933 # batch 205 summa-be32 value
misshru = missgru  # could be different

testing = False
if testing:
    top_fold = '/Users/amedin/Research/USask/test_py/'
    method_name = 'sundials_1en8'
else:
    import multiprocessing as mp
    import sys
    # The first input argument specifies the run where the files are
    method_name = sys.argv[1] # sys.argv values are strings by default so this is fine (sundials_1en8 or be64)

ncdir        = top_fold + 'summa-' + method_name + '_nocat'
file_pattern = 'run1_G*_timestep.nc'
ctdir        = top_fold + 'summa-' + method_name

# get list of split summa output files (hardwired pattern)
outfilelist0 = glob((ncdir+'/'+file_pattern))
outfilelist0.sort()

# -- functions
def get_stat(g,catby_num,outfilelist0,ctdir):
    outfilelist = outfilelist0[(catby_num*g):(catby_num*(g+1))]
    gru_num = 0
    hru_num = 0
    subset0 = outfilelist[0].split('/')[-1].split('_')[1]
    subset1 = outfilelist[-1].split('/')[-1].split('_')[1]
    out_name = 'run1_'+subset0[0:7]+subset1[7:14]+'_timestep.nc' # will fail if GRU numbers are more than 6 digits

    for file in outfilelist:
        f = nc.Dataset(file)
        gru_num = gru_num+len(f.dimensions['gru'])
        hru_num = hru_num+len(f.dimensions['hru'])
        # extract the subset IDs

        # write output
    with nc.Dataset(outfilelist[0]) as src:
        with nc.Dataset(ctdir+'/'+out_name, "w") as dst:
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
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                # Note here the variable dimension name is the same, but size has been updated for gru and hru.

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
            for i,file in enumerate(outfilelist):

                print("combining file %d %s" % (i,file))
                # f = nc.Dataset(os.path.join(ncdir, file))
                f = nc.Dataset(file)
                for j in range(gru_vars_num):
                    gru_var_name = gru_vars[j][0]
                    dim_index = gru_vars[j][1]
                    data=f[gru_var_name][:]
                    if i == 0:
                        Dict[gru_var_name]=data
                    else:
                        Dict[gru_var_name]=np.concatenate((Dict[gru_var_name],data),axis=dim_index)

                for j in range(hru_vars_num):
                    hru_var_name = hru_vars[j][0]
                    dim_index = hru_vars[j][1]
                    data=f[hru_var_name][:]
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
                new_index = np.append(dst["gru"].values,missgru)
                dst.reindex({"gru": new_index})
                dst.sel(gru=missgru)["gruId"] = missgru

                new_index = np.append(dst["hru"].values,misshru)
                dst.reindex({"hru": new_index})
                dst.sel(gru=misshru)["hruId"] = misshru

            # Temporarily create gruId from hruId
            #if gru_num == hru_num:
            #    gruId = dst.createVariable('gruId', dst['hruId'].datatype, ('gru',))
            #    gruId.long_name = "ID of group of response unit (GRU)"
            #    gruId.units = dst['hruId'].units
            #    dst.variables['gruId'][:] = dst.variables['hruId'][:]
            #else:
            #    print('Warning: gruId variable cannot be created since it has different size from hruId')

        print("wrote output: %s" % (ctdir+'/'+out_name))

    return #nothing
# -- end functions


if testing:
    # -- no parallel processing
    for g in range(0,int(len(outfilelist0)/catby_num)):
        get_stat(g,catby_num,outfilelist0,ctdir)
else:
    # -- start parallel processing
    ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
    if __name__ == "__main__":
        pool = mp.Pool(processes=ncpus)
        results = [pool.apply_async(get_stat, args=(g,catby_num,outfilelist0,ctdir)) for g in range(0,int(len(outfilelist0)/catby_num))]
        dojob = [p.get() for p in results]
        pool.close()
    # -- end parallel processing


