"""Check the NCAR/summa develop port of J() against the validated summa-spotpy build.

One J() evaluation of the basin's best calibrated parameters over its calibration run
(settings/SUMMA/fileManager.txt), in one of three variants:
  old       summa-spotpy libsumma + settings/<basin>_config.toml              (reference)
  dev_ts    develop libsumma + settings/<basin>_config_develop_writeTS.toml    (must match old)
  dev_nots  develop libsumma + settings/<basin>_config_develop_noTS.toml       (no time-series output)
Usage: python validate_basin.py <ross|athabasca> <old|dev_ts|dev_nots> [tag]
Prints one RESULT line: basin, variant, KGE, wall time, size of output files written.
"""
import json, os, re, sys, time
from pathlib import Path

ROOT = Path(os.environ.get("WORK_ROOT", "<WORK_ROOT>"))   # folder holding summa-spotpy/ and summa-develop/
SPOT = ROOT / "summa-spotpy"
DEV_LIB = ROOT / "summa-develop" / "build" / "cmake_build" / "libsumma.so"
PORT = SPOT / "my_work_athabasca" / "scripts" / "develop_port"      # holds the SUMMA_LIB-aware interface
basin_name, variant = sys.argv[1], sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else os.environ.get("SLURM_ARRAY_JOB_ID", "local")
basin = SPOT / f"my_work_{basin_name}"
best_run = {"ross": "asyncdds_mpi_21008912", "athabasca": "asyncdds_mpi_21150407"}[basin_name]
cfg_name = {"old": f"{basin_name}_config.toml",
            "dev_ts": f"{basin_name}_config_develop_writeTS.toml",
            "dev_nots": f"{basin_name}_config_develop_noTS.toml"}[variant]

p = json.load(open(basin / "results" / best_run / "summary.json"))["best_params"]
work = basin / "output" / f"devcheck_{tag}" / variant
work.mkdir(parents=True, exist_ok=True)
# develop renamed one output variable (scalarSnowSublimation -> scalarGroundSublimation) and
# stops on the old name, so develop runs read a converted copy of outputControl.txt
RENAMED = {"scalarSnowSublimation": "scalarGroundSublimation"}
oc_dev = basin / "settings" / "SUMMA" / "outputControl_develop.txt"
if variant != "old":
    txt = (basin / "settings" / "SUMMA" / "outputControl.txt").read_text()
    for a, b in RENAMED.items():
        txt = re.sub(rf"^(\s*){a}\b", rf"\g<1>{b}", txt, flags=re.M)
    if not oc_dev.exists() or oc_dev.read_text() != txt:
        oc_dev.write_text(txt)
def fm_line(l):
    k = l.strip()
    if k.startswith("outputPath"):
        return f"outputPath '{work}/'"
    if variant != "old" and k.startswith("outputControlFile"):
        return "outputControlFile 'outputControl_develop.txt'"
    return l
fm = [fm_line(l) for l in (basin / "settings" / "SUMMA" / "fileManager.txt").read_text().splitlines()]
(work / "fileManager.txt").write_text("\n".join(fm) + "\n")
cfg = re.sub(r'^hfabric_path\s*=.*$', f'hfabric_path = "{work}/"',
             (basin / "settings" / cfg_name).read_text(), flags=re.M)
(work / "config.toml").write_text(cfg)
if not (work / "topology.nc").exists():
    os.symlink(basin / "settings" / "mizuRoute" / "topology.nc", work / "topology.nc")
os.environ["SUMMA_MASTER"] = str(work / "fileManager.txt")
os.environ["SUMMA_CONFIG"] = str(work / "config.toml")
if variant == "old":
    sys.path.insert(0, str(SPOT / "calibration"))
else:
    os.environ["SUMMA_LIB"] = str(DEV_LIB)
    sys.path.insert(0, str(PORT))
os.chdir(work)
from summa_python_interface import J, PARAM_NAMES, LIB_PATH

t0 = time.time()
kge = J([p[n] for n in PARAM_NAMES])
wall = time.time() - t0
out_mb = sum(f.stat().st_size for f in work.glob("*.nc") if f.name != "topology.nc") / 1e6
print(f"RESULT {basin_name:9s} {variant:9s} KGE = {kge:.4f}   wall {wall/60:.2f} min   "
      f"output files {out_mb:.1f} MB   lib {LIB_PATH}", flush=True)
