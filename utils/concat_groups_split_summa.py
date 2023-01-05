# concatenate the outputs of a split domain summa run into fewer groups
# written originally by Manab Saharia, updated by Hongli Liu and Andy Wood, and W. Knoben
# modified by A. Van Beusekom (2023)
# Usage: python concat_split_summa.py

import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
from glob import glob
import netCDF4 as nc
import numpy as np

method_name = 'sundials_1en8'
catby_num    = 2 #number of files to cat into one
top_fold     = '/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/'
#top_fold     = '/Users/amedin/Research/USask/test_py/'

ncdir        = top_fold + 'summa-' + method_name
file_pattern = 'run1_G*_timestep.nc'
ctdir        = top_fold + 'summa-' + method_name + '_cat'

# get list of split summa output files (hardwired pattern)
outfilelist0 = glob((ncdir+'/'+file_pattern))
outfilelist0.sort()

# make new outlist of catby_num
for i in range(0,int(len(outfilelist0)/catby_num)):
	outfilelist = outfilelist0[catby_num*i:catby_num*(i+1)]
	print((i,len(outfilelist0)/catby_num,outfilelist))
	# count the number of gru and hru
	gru_num = 0
	hru_num = 0
	subset0 = outfilelist[0].split('/')[-1].split('_')[1]
	subset1 = outfilelist[-1].split('/')[-1].split('_')[1]
	out_name = 'run1_'+subset0[0:7]+subset1[7:14]+'_timestep.nc' # will fail if GRU numbers are more than 6 digits

	for file in outfilelist:
		# f = nc.Dataset(os.path.join(ncdir, file))
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

			# Temporarily create gruId from hruId
			#if gru_num == hru_num:
			#    gruId = dst.createVariable('gruId', dst['hruId'].datatype, ('gru',))
			#    gruId.long_name = "ID of group of response unit (GRU)"
			#    gruId.units = dst['hruId'].units
			#    dst.variables['gruId'][:] = dst.variables['hruId'][:]
			#else:
			#    print('Warning: gruId variable cannot be created since it has different size from hruId')

	print("wrote output: %s" % (ctdir+'/'+out_name))

print('Done')
