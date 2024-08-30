#!/usr/bin/env python

# python check_bit_4_bit_withTol.py [bench_dir] [test_dir] [tol]

import sys
import numpy as np
import xarray as xr
from pathlib import Path

PRINT_MESSAGE = 'Test {0: >3} - Filename: {1: <64} - Variable: {2: <64} - Mean and Max difference: {3: <8}, {4: <8}'
FINAL_MESSAGE = 'Variable: {0: >20} - HRU: {1: <8} - Max: {2: <25} - Time: {3: <30}'
ERROR_MESSAGE = 'Found differences in test cases: {}'
SUCCESS_MESSAGE = 'SUCCESSFUL TEST RUN: All files match within tolerence!'
USAGE = 'USAGE: python check_bit_4_bit_withTol.py BENCHMARK_DIRECTORY TEST_DIRECTORY TOLERENCE'

#if len(sys.argv) != 4:
#    print(USAGE)
#    exit(1)

testing = False
if testing:
    bench_dir = 'summa_be' #sys.argv[1]
    test_dir = 'summa-sundials_be' #sys.argv[2]
    tol = 0.1 #sys.argv[3]
else:
    bench_dir = sys.argv[1]
    test_dir = sys.argv[2]
    tol = sys.argv[3]

bench_files = sorted(list(Path(bench_dir).glob('**/*.nc')))
test_files = sorted(list(Path(test_dir).glob('**/*.nc')))


assert len(bench_files) == len(test_files), \
    'Found {} files but need {}!'.format(len(test_files), len(bench_files))


def rem(li1, li2):
    return list(set(li1) - set(li2))


def compute_diffs(bench_files, test_files):
    all_diffs_mean = []
    all_times_max = []
    all_diffs_max = []
    all_tots_mean = []
    all_tots_max  = []
    for i, (f1, f2) in enumerate(zip(bench_files, test_files)):
        ds1, ds2 = xr.open_dataset(str(f1)), xr.open_dataset(str(f2))
        #m = (ds1['scalarSWE'].sel(hru = 71003766)- ds2['scalarSWE'].sel(hru = 71003766))
        m = (ds1['scalarSWE'].sel(time = '1981-05-02T14:00:00.000013408')- ds2['scalarSWE'].sel(time = '1981-05-02T14:00:00.000013408'))
        m.plot()
        diff = (np.fabs(ds1 - ds2))
        # get rid of gru dimension, assuming they are same as the often are
        diff = diff.drop_vars(['hruId','gruId','wallClockTime'])
        m = diff.drop_dims('hru')
        m = m.rename({'gru': 'hru'})
        diff = diff.drop_dims('gru')
        diff = xr.merge([diff,m])
        # take means and maxes
        diff_mean = diff.mean(dim='time')
        time_max = diff.idxmax(dim='time')
        #for v in rem(list(diff.variables.keys()), list(diff.dims.keys())):
        #    m = diff[v].sel(hru = np.int(diff[v].idxmax(dim="hru").values[0]))
        #    print(v, m.idxmax().values,np.int(diff[v].idxmax(dim="hru").values[0]),diff[v].max(dim="hru").values[0])
        diff_max = diff.max(dim='time')
        tot_mean = diff_mean.mean(dim='hru')
        tot_max = diff_max.max(dim='hru')
        #for v in rem(list(diff_mean.variables.keys()), list(diff_mean.dims.keys())):
            #print(PRINT_MESSAGE.format(i, str(f1).split('/')[-1], v, tot_mean[v].values, tot_max[v].values))
        all_diffs_mean.append(diff_mean)
        all_times_max.append(time_max)
        all_diffs_max.append(diff_max)
        all_tots_mean.append(tot_mean)
        all_tots_max.append(tot_max)
    return all_diffs_mean, all_times_max, all_diffs_max, all_tots_mean, all_tots_max,i


all_diffs_mean, all_times_max, all_diffs_max, all_tots_mean, all_tots_max,i = compute_diffs(bench_files, test_files)

combined_mean = xr.concat(all_diffs_mean, dim='hru')
combined_mean.to_netcdf("all_diffs_mean.nc")

combined_time = xr.concat(all_times_max, dim='hru')
combined_max = xr.concat(all_diffs_max, dim='hru')
combined_max.to_netcdf("all_diffs_max.nc")
max_ind = combined_max.idxmax(dim='hru')

combined_tmean = xr.concat(all_tots_mean, dim='hru')
combined_tmax = xr.concat(all_tots_max, dim='hru')

#assert np.fabs(combined_tmean).max(dim='hru') - tol*(i-1) <= 0.0, \
#    ERROR_MESSAGE.format(np.argwhere(np.asarray(combined_tmean)-tol*(i-1)> 0.0))

#print(SUCCESS_MESSAGE)
print(combined_tmean)


# HRU - 71000000 is the GRU if a 71* HRU, gruId, hruId is the first HRU, if value is first HRU then all values same
for v in rem(list(combined_max.variables.keys()), list(combined_time.dims.keys())):
    print(FINAL_MESSAGE.format(v, np.int(max_ind[v].values), combined_max[v].sel(hru = np.int(max_ind[v].values)).values, combined_time[v].sel(hru = np.int(max_ind[v].values)).values))

if testing:
    bench_files = "run_3766_summa_be.nc"
    f1 = bench_files

    #test_files = ["run_3766_summa_be_conv.nc","run_3766_summa-sundials_beclosed_new.nc","run_3766_summa-sundials_closed_new.nc","run_3766_summa-sundials_closed.nc","run_3766_summa-sundials_enthal_new.nc","run_3766_summa-sundials_enthal_new_noderiv.nc","run_3766_summa-sundials_enthal.nc","run_3766_summa-sundials_beclosed_new_noCp.nc","run_3766_summa-sundials_closed_new_noCp.nc"]
    #test_files = ["run_3766_summa-sundials_beclosed_new.nc","run_3766_summa-sundials_beclosed_new_noCp.nc","run_3766_summa-sundials_beclosed_new_noCp_oldJac.nc"]
    test_files = ["run_3766_sun_updateCp_withDeriv.nc","run_3766_be_updateCp_withDeriv.nc","run_3766_be_original_withDeriv.nc","run_3766_be_orginal.nc"]
    vars = ['averageRoutedRunoff','scalarTotalET','scalarTotalSoilWat','scalarSWE','scalarCanopyWat']

    for j, v in enumerate(vars):
        for i, f2 in enumerate(test_files):
            ds1, ds2 = xr.open_dataset(f1), xr.open_dataset(f2)
            m = ds1[v]- ds2[v]
            m.plot()

        plt.gca().legend(test_files)
        plt.show()

    test_files = ["run_3766_nrgWatTermOrderSwitch_commita2c16c8c.nc","run_3766_nrgWatTermOrderSame_commitec1f42b4.nc"]
    vars = ['averageRoutedRunoff','scalarTotalET','scalarTotalSoilWat','scalarSWE','scalarCanopyWat']

    for j, v in enumerate(vars):
        for i, f2 in enumerate(test_files):
            ds1, ds2 = xr.open_dataset(f1), xr.open_dataset(f2)
            m = ds1[v]- ds2[v]
            m.plot()

        plt.gca().legend(test_files)
        plt.show()

