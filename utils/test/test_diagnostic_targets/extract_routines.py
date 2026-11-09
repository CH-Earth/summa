#!/usr/bin/env python3
"""Compile the actual scalarizer/checker in a light harness, not a SUMMA build."""
from pathlib import Path
import re
import sys
root, out = map(Path, sys.argv[1:])

def routine(path: str, name: str, kind: str) -> str:
    source = (root / path).read_text()
    pattern = rf'^  (?:pure )?{kind} {name}\b.*?^  end {kind} {name}[^\n]*'
    matches = re.findall(pattern, source, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if len(matches) != 1:
        raise RuntimeError(f'Expected exactly one {name} in {path}, found {len(matches)}')
    return matches[0] + '\n'

simulation = 'build/source/driver/summa_simulation.f90'
header = '''module summa_type
 use calibration_types, only: calib_info
 type :: config_info
   type(calib_info) :: calib
   logical :: write_timeseries=.false.
 end type
end module
module summa_simulation
 use nr_type, only: i4b,rkind,lgt
 use summa_type, only: config_info
 contains
'''
body = routine(simulation, 'n_calibration_targets', 'function')
body += routine(simulation, 'scalarize_objectives', 'function')
(out/'scalarizer_under_test.f90').write_text(header+body+'end module\n')
body = routine('build/source/calibration/summa_parameter_sampling.f90', 'check_search_settings', 'subroutine')
(out/'settings_under_test.f90').write_text('''module settings_under_test
 use nr_type, only: i4b,rkind,lgt
 use summa_type, only: config_info
 contains
'''+body+'end module\n')
body = routine(simulation, 'score_target', 'subroutine')
(out/'scorer_under_test.f90').write_text('''module scorer_under_test
 use nr_type, only: i4b,rkind,lgt
 use summa_type, only: config_info
 use globalData, only: integerMissing,realMissing,iulog
 contains
'''+body+'end module\n')
