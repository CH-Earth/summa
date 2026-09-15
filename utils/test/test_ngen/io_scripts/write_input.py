#!/usr/bin/env python3
"""
Write input files for NGEN from SUMMA fileManager and NGEN forcing files.
Usage:
  python3 write_input.py summa_fileManager.nc ngen_forcing.nc restart_freq output_folder 
"""
import argparse
import xarray as xr
import numpy as np
import re
import os

def _parse_id(fid):
    # try direct int, else extract trailing digits, else first digits, else return original string
    try:
        return int(fid)
    except Exception:
        s = str(fid).strip()
        m = re.search(r'(\d+)$', s) or re.search(r'(\d+)', s)
        return int(m.group(1)) if m else s
    
def main():
    p = argparse.ArgumentParser()
    p.add_argument("fileManager", help="summa_fileManager.nc (SUMMA fileManager file)")
    p.add_argument("ngen_forcing", help="ngen_forcing.nc (contains target HRU ids/order)")
    p.add_argument("restart_freq", help="restart frequency string, y for every year, m for month, d for day, e for end, n for never")
    p.add_argument("output_folder", help="output folder for NGEN input files")
    args = p.parse_args() 

    # forcing ids (use as labels for the new hru coordinate)
    f_ds = xr.open_dataset(args.ngen_forcing)
    f_var = "ids"
    if f_var not in f_ds:
        raise RuntimeError(f"forcing file missing '{f_var}' variable")
    f_ids_raw = np.array([str(x) for x in f_ds[f_var].values])
    f_ids = np.array([_parse_id(fid) for fid in f_ids_raw])

    # ensure output folder exists
    os.makedirs(args.output_folder, exist_ok=True)
    current_dir = os.getcwd()
    # get directory starting from "extern"
    idx = current_dir.find("extern")
    extern_path = current_dir[idx:] if idx != -1 else current_dir

    # write a simple &parameters namelist block for each HRU id (stub)
    for i, id_ in enumerate(f_ids):
        out_path = os.path.join(args.output_folder, f"cat-{id_}.input")
        with open(out_path, "w") as fh:
            fh.write("&parameters\n")
            fh.write(f'  file_manager          = "./{extern_path}/{args.fileManager}"\n')
            fh.write(f'  attrib_file_HRU_order = {str(i+1)}\n')
            fh.write(f'  restart_print_freq    = "{args.restart_freq}"\n')
            fh.write("/\n")
        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()