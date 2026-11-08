#!/usr/bin/env python3
"""Extract the unchanged production calibration-file definition for metadata tests."""
from pathlib import Path
import re
import sys
root, out = map(Path, sys.argv[1:])
src = (root/'build/source/calibration/calibration_output.f90').read_text()
def routine(name, kind):
    found = re.findall(rf'^  {kind} {name}\b.*?^  end {kind} {name}[^\n]*', src, re.M|re.S|re.I)
    if len(found) != 1:
        raise RuntimeError(f'Expected exactly one {name}, found {len(found)}')
    return found[0]+'\n'
(out/'calibration_output_under_test.f90').write_text('''module calibration_output_under_test
 use nr_type, only:i4b,rkind,lgt
 use netcdf
 use calibration_types, only:target_info
 use parameter_search, only:parameter_spec,param_physical,param_spatial_multiplier,param_ratio
 implicit none
 contains
''' + routine('create_calibration_output','subroutine') + routine('joined_field','function') + 'end module\n')
