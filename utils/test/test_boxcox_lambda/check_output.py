#!/usr/bin/env python3
"""Check Box-Cox metadata in a real SUMMA calibration NetCDF file.

Optional; requires numpy and netCDF4. Example:
  python3 check_output.py calibration.nc --expect streamflow_nse_boxcox=0.2
This checks saved metadata, not the numerical metric values or the underlying flow.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('calibration_file', type=Path)
    parser.add_argument('--expect', action='append', default=[], metavar='TARGET=LAMBDA',
                        help='Check a target-specific value; may be repeated.')
    args = parser.parse_args()
    try:
        import numpy as np
        from netCDF4 import Dataset, chartostring
    except ImportError as exc:
        print(f'This optional check requires numpy and netCDF4: {exc}', file=sys.stderr)
        return 2
    try:
        expected = {}
        for text in args.expect:
            name, value = text.rsplit('=', 1)
            value = float(value)
            if not name or not np.isfinite(value):
                raise ValueError(f'Invalid expected exponent: {text}')
            expected[name] = value
        with Dataset(args.calibration_file) as ds:
            name_var = ds.variables['target_name']
            raw = name_var[:]
            if raw.dtype.kind == 'S' and raw.dtype.itemsize == 1:
                raw = np.moveaxis(raw, name_var.dimensions.index('target'), 0)
                raw = chartostring(raw)
            names = [x.decode().strip(' \0') if isinstance(x, bytes) else str(x).strip(' \0') for x in raw]
            transforms = [x.strip() for x in ds.variables['objective'].getncattr('obs_transform').split(',')]
            var = ds.variables['target_boxcox_lambda']
            if var.dimensions != ('target',):
                raise ValueError(f'Unexpected lambda dimensions: {var.dimensions}')
            values = np.ma.asarray(var[:], dtype=float)
            if len(names) != len(transforms) or values.shape != (len(names),):
                raise ValueError('Target names, transformations and lambda dimensions differ')
            missing = np.ma.getmaskarray(values)
            for j, (name, transform) in enumerate(zip(names, transforms)):
                if transform == 'boxcox':
                    if missing[j] or not np.isfinite(values[j]):
                        raise ValueError(f'Missing/nonfinite exponent for Box-Cox target {name}')
                    print(f'  {name}: boxcox_lambda = {float(values[j]):.17g}')
                elif not missing[j]:
                    raise ValueError(f'Expected missing lambda for non-Box-Cox target {name}')
            for name, value in expected.items():
                j = names.index(name)
                if transforms[j] != 'boxcox' or missing[j] or not np.isclose(values[j], value, rtol=1e-13, atol=0):
                    raise ValueError(f'Unexpected exponent for {name}; expected {value}')
        print(f'PASS: Box-Cox metadata in {args.calibration_file}')
        return 0
    except (OSError, KeyError, ValueError, AttributeError) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
