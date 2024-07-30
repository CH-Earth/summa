#!/usr/bin/env python
''' Create a vector cold state file for SUMMA from constant values'''
#
# Author:    Andy Wood, Feb 2017
#
# Notes:  quick n dirty to generate constant initial states across a domain
#         all values hardwired, just gets HRU index from an existing parameter file
#         improvements:  could read another cold state file to get list of variables
#                        to populate; or read a metadata dictionary (names, types, etc)
#
#         no mapping required here -- but one could map another resolution vals, similar
#         to the param mapping scripts
#
#         check: nSoil and nSnow might have to be 'int' instead of 'int64' in output
#
# Requirements:  run with a python (eg miniconda) 2.7 that includes netCDF4
#
# =========================================================================

import sys
import os
import time
import getopt
import numpy as np
import netCDF4 as nc4
#import xarray as xr

########################################################################
#                                Subroutines                           #
########################################################################

def getNetCDFData(fn, varname):
    """Read <varname> variables available to be mapped from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
#    ds = xr.open_dataset(fn)
#    data = ds[varname]
    return data

def getOutputPolyIDs(nc_file):
    outPolyIDs  = getNetCDFData(nc_file, 'hruId')    
    print("read output outPolyIds ('hruId') from example domain file")
    return outPolyIDs


# write variables to netcdf output file
def writeNC_state_vars(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim, 'hru',),fill_value='-999.0')    
    ncvar[:] = newVarVals   # store data in netcdf file


# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn,  scalarv, midSoil, midToto, ifcToto, hrus, hru_type):    
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("writing output file")
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_hru = nc_out.createDimension('hru', len(hrus))
    dim_scalarv = nc_out.createDimension('scalarv', scalarv)
    dim_midSoil = nc_out.createDimension('midSoil', midSoil)
    dim_midToto = nc_out.createDimension('midToto', midToto)
    dim_ifcToto = nc_out.createDimension('ifcToto', ifcToto)    

    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20  # EC
        dim_str = nc_out.createDimension('strlen', max_strlen)
        hruId = nc_out.createVariable('hruId', 'S1', ('hru', 'strlen'),fill_value='-999')  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))         
    elif hru_type == 'int64':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'i8', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')

    elif hru_type == 'int':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'int', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')

    else:
        # not recognized
        sys.exit("ERROR, hru_type not recognized: must be str, int64, or int")

    # add attribute    
    hruId.long_name = 'USGS HUC12 ID'

    return nc_out
    # leave netcdf file open


############################################
#                  Main                    #
############################################
use = '''
Usage: %s -[h] <existing_inputfile_with_hruId> <output_netCDF> <hru_type (int|int64|str)>
        -h  help
'''
if __name__ == '__main__':

    def usage():
        sys.stderr.write(use % sys.argv[0])
        sys.exit(1)
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], 'h')
    except getopt.error:
        usage()

    verbose   = False
    grid_info = False
    proj_info = True
    for (opt,val) in opts:
        if opt == '-h':
            usage()
        elif opt == '-v':
            verbose = True
        else:
            raise OptionError(opt)
            usage()

    if len(args) == 3:
        nc_example_name = args[0]   # template file (other params)
        nc_out_name = args[1]   # output cold-state file
        hru_type = args[2]   # 'int' or 'string'
        # nc_in_name  = args[4]   # existing cold-state file
 
        # hardwired to forcing formats (hru index rather than grid)
        outPolyIDs=getOutputPolyIDs(nc_example_name)        
        nOutPolygons = len(outPolyIDs)

        # === now start to create the cold state variables using the variable template ===

        # settings (hardwire for now - later read from a config file -- 4 layer
        #scalarv = 1
        #midSoil = 4
        #midToto = 4
        #ifcToto = 5
        #dT = 10800 # timestep of forcings in seconds
        #lyrDepth  = [0.1, 0.3, 0.6, 1.0]
        #lyrHeight = [0.0, 0.1, 0.4, 1.0, 2.0]

        # settings 3 layer
        scalarv = 1
        midSoil = 3
        midToto = 3
        ifcToto = 4
        dT = 10800 # timestep of forcings in seconds
        #dT = 3600 # timestep of forcings in seconds
        lyrDepth  = [0.1, 0.3, 0.6]
        lyrHeight = [0.0, 0.1, 0.4, 1.0]

        # initialize netcdf file by storing dimensions and hru variable
        nc_out = writeNC_dims(nc_out_name, scalarv, midSoil, midToto, ifcToto,
                              outPolyIDs, hru_type)

        # === now loop through variables and write
        #  this could be done by looping through the input state file and xferring values
        #  could also read vars and 

        # nSoil, nSnow
        newVarVals = np.full((1,nOutPolygons), midSoil, dtype='f8')
        writeNC_state_vars(nc_out, 'nSoil', 'scalarv', 'f8', newVarVals)
        newVarVals = np.zeros((1,nOutPolygons), dtype='f8')
        writeNC_state_vars(nc_out, 'nSnow', 'scalarv', 'f8', newVarVals)                

        # dT
        newVarVals = np.full((1,nOutPolygons), dT)        
        writeNC_state_vars(nc_out, 'dt_init', 'scalarv', 'f8', newVarVals)

        # SWE, SnowDepth, SfcMeltPond, SnowAlbedo, CanopyLiq, CanopyIce
        newVarVals = np.zeros((1,nOutPolygons))
        writeNC_state_vars(nc_out, 'scalarSWE', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarSnowDepth', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarSfcMeltPond', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarSnowAlbedo', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarCanopyLiq', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarCanopyIce', 'scalarv', 'f8', newVarVals)        

        # CanairTemp, CanopyTemp
        newVarVals = np.full((1,nOutPolygons), 283.16)        
        writeNC_state_vars(nc_out, 'scalarCanairTemp', 'scalarv', 'f8', newVarVals)
        writeNC_state_vars(nc_out, 'scalarCanopyTemp', 'scalarv', 'f8', newVarVals)

        # AquiferStorage
        newVarVals = np.full((1,nOutPolygons), 1.0)        
        writeNC_state_vars(nc_out, 'scalarAquiferStorage', 'scalarv', 'f8', newVarVals)

        # layer MatricHead
        newVarVals = np.full((midSoil,nOutPolygons), -1.0)
        writeNC_state_vars(nc_out, 'mLayerMatricHead', 'midSoil', 'f8', newVarVals)

        # layer Temp
        newVarVals = np.full((midToto,nOutPolygons), 283.16)
        writeNC_state_vars(nc_out, 'mLayerTemp', 'midToto', 'f8', newVarVals)

        # layer VolFracLiq
        newVarVals = np.full((midToto,nOutPolygons), 0.2)
        writeNC_state_vars(nc_out, 'mLayerVolFracLiq', 'midToto', 'f8', newVarVals)

        # layer VolFracIce
        newVarVals = np.full((midToto,nOutPolygons), 0.0)
        writeNC_state_vars(nc_out, 'mLayerVolFracIce', 'midToto', 'f8', newVarVals)

        # layer Depth, Height
        newVarVals = np.full((nOutPolygons,midToto), lyrDepth).transpose()
        writeNC_state_vars(nc_out, 'mLayerDepth', 'midToto', 'f8', newVarVals)
        newVarVals = np.full((nOutPolygons,ifcToto), lyrHeight).transpose()
        writeNC_state_vars(nc_out, 'iLayerHeight', 'ifcToto', 'f8', newVarVals)        

        nc_out.close()

    else:
        usage()


        # code for reading input data variable & attributes
        # to pass to processing and write routines (instead of hardwired calls above
#        f = nc4.Dataset(nc_infl, 'r')
#        var_in = f.variables[varname]
#        attNames = []
#        attContents = []
#        attr = var_in.ncattrs()  # get attributes
#        for n in range(0,len(attr)):
#            attNames.extend([attr[n]])
#            attContents.extend([var_in.getncattr(attr[n])])
