#!/usr/bin/env python
import sys
import numpy as np
import xarray as xr
from pathlib import Path

PRINT_MESSAGE = 'Test {0: >3} - Filename: {1: <64} - Total difference: {2: <8}'
ERROR_MESSAGE = 'Found differences in test cases: {}'
SUCCESS_MESSAGE = 'SUCCESSFUL TEST RUN: All files match!'
USAGE = 'USAGE: ./check_bit_4_bit.py BENCHMARK_DIRECTORY TEST_DIRECTORY'

if len(sys.argv) != 3:
    print(USAGE)
    exit(1)

bench_dir = sys.argv[1]
test_dir = sys.argv[2]

bench_files = sorted(list(Path(bench_dir).glob('**/*.nc')))
test_files = sorted(list(Path(test_dir).glob('**/*.nc')))


assert len(bench_files) == len(test_files), \
    'Found {} files but need {}!'.format(len(test_files), len(bench_files))


def rem(li1, li2):
    return list(set(li1) - set(li2))


def compute_diffs(bench_files, test_files):
    all_diffs = []
    all_tots = []
    for i, (f1, f2) in enumerate(zip(bench_files, test_files)):
        ds1, ds2 = xr.open_dataset(str(f1)), xr.open_dataset(str(f2))
        diff = (ds1 - ds2).sum(dim='time')
        tot = 0.0
        for v in rem(list(diff.variables.keys()), list(diff.dims.keys())):
            tot += np.sum(diff[v].values)
        all_diffs.append(diff)
        all_tots.append(tot)
        print(PRINT_MESSAGE.format(i, str(f1).split('/')[-1], tot))
    return all_diffs, all_tots


all_diffs, all_tots = compute_diffs(bench_files, test_files)

assert np.sum(np.absolute(all_tots)) == 0.0, \
    ERROR_MESSAGE.format(np.argwhere(np.asarray(all_tots) != 0))

print(SUCCESS_MESSAGE)
