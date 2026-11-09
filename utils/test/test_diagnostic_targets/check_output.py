#!/usr/bin/env python3
"""Inspect a REAL SUMMA calibration result after adding the example targets.

Requires numpy and netCDF4. This is an output-schema/validity check, not a
recalculation of the metrics from the underlying model time series.
"""
import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    from netCDF4 import Dataset, chartostring
except ImportError as exc:
    raise SystemExit('This optional output check needs numpy and netCDF4 in your Python environment.') from exc


def check_output(path: Path, active_name: str, diagnostic_name: str) -> None:
    with Dataset(path) as ds:
        names_var = ds.variables['target_name']
        if names_var.dtype.kind == 'S' and names_var.dtype.itemsize == 1:
            names = [str(x).strip() for x in chartostring(names_var[:])]
        else:
            names = [str(x).strip() for x in names_var[:]]
        weights = np.asarray(ds.variables['target_weight'][:], dtype=float)
        if weights.shape != (len(names),):
            raise ValueError('target_weight and target_name lengths differ')
        if not np.isfinite(weights).all() or (weights < 0).any() or not (weights > 0).any():
            raise ValueError('invalid optimization weights')
        i_active, i_diag = names.index(active_name), names.index(diagnostic_name)
        if weights[i_active] != 1.0 or weights[i_diag] != 0.0:
            raise ValueError('expected KGE weight 1.0 and diagnostic weight 0.0')
        var = ds.variables['objective']
        if set(var.dimensions) != {'sample', 'target'}:
            raise ValueError(f'unexpected objective dimensions: {var.dimensions}')
        obj = np.ma.asarray(var[:], dtype=float)
        obj = np.moveaxis(obj, var.dimensions.index('sample'), 0)
        values = obj.filled(np.nan)
        sample = np.ma.asarray(ds.variables['sample'][:])
        completed = (~np.ma.getmaskarray(sample)) & (sample.filled(-1) > 0)
        failed = np.ma.asarray(ds.variables['trial_failed'][:]).filled(1) != 0
        success = completed & ~failed
        if not success.any():
            raise ValueError('no successful completed trials to inspect')
        if not np.isfinite(values[success][:, weights > 0]).all():
            raise ValueError('successful trial contains a missing active objective')
        diagnostic_missing = int(np.sum(~np.isfinite(values[success, i_diag])))
        print(f'PASS: {path}')
        print(f'  targets: {dict(zip(names, weights.tolist()))}')
        print(f'  completed: {completed.sum()}, successful: {success.sum()}, failed: {(completed & failed).sum()}')
        print(f'  undefined diagnostics in successful trials: {diagnostic_missing}')
        print('  Diagnostic missing values do not mark the whole trial as failed.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('calibration_file', type=Path)
    parser.add_argument('--active', default='streamflow_kge')
    parser.add_argument('--diagnostic', default='streamflow_nse_log')
    args = parser.parse_args()
    try:
        check_output(args.calibration_file, args.active, args.diagnostic)
    except (OSError, KeyError, ValueError) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 1
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
