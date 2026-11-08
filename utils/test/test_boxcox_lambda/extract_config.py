#!/usr/bin/env python3
"""Extract unchanged production parsers/legacy adapter into lightweight test modules."""
from pathlib import Path
import re
import sys
root, out = map(Path, sys.argv[1:])
def routine(path, name):
    text = (root / path).read_text()
    found = re.findall(rf'^  subroutine {name}\b.*?^  end subroutine {name}[^\n]*',
                       text, re.M | re.S | re.I)
    if len(found) != 1:
        raise RuntimeError(f'Expected exactly one {name}, found {len(found)}')
    return found[0] + '\n'
path = 'build/source/calibration/calibration_config.f90'
body = ''.join(routine(path, name) for name in [
    'parse_calibration_settings', 'parse_calibration_targets', 'parse_target_spatial_unit'])
(out / 'parsers_under_test.f90').write_text('''module parsers_under_test
 use nr_type, only:i4b,i8b,rkind,lgt
 use calibration_types, only:calib_info,target_info
 use tomlf_all, only:toml_table,toml_key,get_value
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 implicit none
 contains
''' + body + 'end module\n')
body = routine('build/source/driver/summa_config.f90', 'set_default_calibration_target')
(out / 'legacy_under_test.f90').write_text('''module legacy_under_test
 use nr_type, only:i4b,rkind
 use calibration_types, only:calib_info
 implicit none
 type :: obs_info
   character(:),allocatable :: obs_path,obs_file,vname_obsflow
 end type
 type :: config_info
   type(calib_info) :: calib
   type(obs_info) :: obs
 end type
 contains
''' + body + 'end module\n')
