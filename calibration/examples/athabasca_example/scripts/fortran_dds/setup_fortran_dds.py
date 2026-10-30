"""Create the input files for a run of SUMMA develop's own (Fortran) async DDS on Athabasca, with the
same model setup as the J() runs. Writes everything into my_work_athabasca/fortran_dds/<tag>/:
  fileManager.txt  : copy of settings/SUMMA/fileManager.txt; simStartTime -> --sim-start (the Fortran
                     driver first spins up for one year before it), outputPath -> <tag>/output/, the develop
                     output list (outputControl_develop.txt, as the J() develop runs use), and the
                     localParamInfo copy below
  localParamInfo_fortran.txt (in settings/SUMMA/): copy with theta_sat's lower bound raised to
                     --theta-sat-low (the Fortran DDS aborts the whole job on the first failed SUMMA run,
                     and theta_sat below ~0.383 fails); all other bounds unchanged
  config.toml      : settings/athabasca_config_develop_noTS.toml + the calibration settings:
                     15 parameters, n_samples, log transforms and ordered constraints as in
                     scripts/athabasca_dds_v2_settings.toml. seg_outlet is removed: in subset mode every
                     model run writes a new topology file into hfabric_path, which 49 ranks would share;
                     without it the reach with the largest area (the outlet, 222) is scored.
Usage (from summa-spotpy/my_work_athabasca):
  python scripts/fortran_dds/setup_fortran_dds.py <tag> <n_samples> [--theta-sat-low 0.39]"""
import argparse, csv, re, sys
from pathlib import Path
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

basin = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("tag"); ap.add_argument("n_samples", type=int)
ap.add_argument("--theta-sat-low", type=float, default=0.39)
ap.add_argument("--sim-start", default="2015-09-01 01:00",
                help="simStartTime for the Fortran run; summa_driver_opt.f90 always spins up for the year "
                     "before it with the default parameters, so forcing must start one year earlier")
a = ap.parse_args()
run = basin / "fortran_dds" / a.tag
(run / "output").mkdir(parents=True, exist_ok=True)
(run / "state").mkdir(parents=True, exist_ok=True)
sset = basin / "settings" / "SUMMA"

# ---- localParamInfo copy with the theta_sat lower bound raised ----
lpi = (sset / "localParamInfo.txt").read_text().splitlines()
done = 0
for k, line in enumerate(lpi):
    parts = line.split("|")
    if len(parts) == 4 and parts[0].strip() == "theta_sat":
        lo = parts[2]
        lpi[k] = "|".join([parts[0], parts[1], f"{a.theta_sat_low:{len(lo) - 1}.4f} ", parts[3]])
        done += 1
assert done == 1, "theta_sat line not found exactly once in localParamInfo.txt"
lpi_name = "localParamInfo_fortran.txt"
(sset / lpi_name).write_text("\n".join(lpi) + "\n")

# ---- fileManager copy ----
fm, seen = [], set()
for line in (sset / "fileManager.txt").read_text().splitlines():
    key = line.split()[0] if line.split() else ""
    if key == "settingsPath":
        fm.append(line); line = f"statePath    '{run / 'state'}/'"; seen.add(key)
    elif key == "statePath":
        continue
    elif key == "outputPath":
        line = f"outputPath '{run / 'output'}/'"; seen.add(key)
    elif key == "outputControlFile":
        line = "outputControlFile    'outputControl_develop.txt'"; seen.add(key)
    elif key == "simStartTime":
        line = f"simStartTime         '{a.sim_start}'"; seen.add(key)
    elif key == "globalHruParamFile":
        line = f"globalHruParamFile   '{lpi_name}'"; seen.add(key)
    fm.append(line)
assert seen == {"settingsPath", "simStartTime", "outputPath", "outputControlFile", "globalHruParamFile"}, f"fileManager keys missing: {seen}"
assert (sset / "outputControl_develop.txt").exists(), "settings/SUMMA/outputControl_develop.txt not found"
(run / "fileManager.txt").write_text("\n".join(fm) + "\n")
# statePath: the spinup restart is written to, and every trial reads it from, the same folder
icond = next(l.split("'")[1] for l in fm if l.split() and l.split()[0] == "initConditionFile")
link = run / "state" / icond
if not link.exists():
    link.symlink_to(sset / icond)

# ---- calibration TOML ----
names = [r["name"] for r in csv.DictReader(open(basin / "scripts" / "athabasca_params.csv"))]
v2 = tomllib.loads((basin / "scripts" / "athabasca_dds_v2_settings.toml").read_text())
cfg = (basin / "settings" / "athabasca_config_develop_noTS.toml").read_text()
cfg, nseg = re.subn(r"^seg_outlet\s*=.*$", "# seg_outlet removed for the Fortran DDS (see setup_fortran_dds.py); "
                    "the largest-area reach (222) is scored", cfg, flags=re.M)
assert nseg == 1, "seg_outlet line not found exactly once"
m = re.search(r"^end_date\s*=.*$", cfg, flags=re.M)
assert m and cfg.rfind("[calibration]", 0, m.start()) > cfg.rfind("\n[", 0, m.start()) - 1, "[calibration] end_date not found"
plist = "param_list = [\n" + ",\n".join(f'    "{n}"' for n in names) + "\n]\n"
cfg = cfg[:m.end()] + f"\nn_samples     = {a.n_samples}\nwrite_aligned = false\n" + plist + cfg[m.end():]
cfg += "\n\n[calibration.parameter_transformations]\n"
cfg += "".join(f'{n} = "{t}"\n' for n, t in v2.get("transformations", {}).items())
for oc in v2.get("ordered", []):
    fixed_only = [q for q in oc["parameters"] if q not in names]
    if fixed_only:
        # summa_driver_opt overrides non-calibrated constraint parameters with their localParamInfo default in
        # every HRU, replacing the per-HRU soil-table values; for critSoilTranspire this fails SUMMA's paramCheck
        # (critSoilTranspire < critSoilWilting). J() never overrides them, so such constraints are left out.
        cfg += f"\n# ordered constraint {oc['parameters']} left out: it would override {fixed_only} in every HRU\n"
        continue
    cfg += "\n[[parameter_dependencies.ordered]]\n"
    cfg += "parameters   = [" + ", ".join(f'"{p}"' for p in oc["parameters"]) + "]\n"
    cfg += f"gap_fraction = {oc.get('gap_fraction', 0.0)}\n"
parsed = tomllib.loads(cfg)                       # must be valid TOML
assert parsed["calibration"]["n_samples"] == a.n_samples and len(parsed["calibration"]["param_list"]) == 15
(run / "config.toml").write_text(cfg)

print(f"Fortran DDS run folder: {run}")
print(f"  n_samples {a.n_samples}, 15 parameters, theta_sat lower bound {a.theta_sat_low}, simStartTime {a.sim_start} "
      f"(one-year default-parameter spinup before it)")
print(f"  log: {list(parsed['calibration']['parameter_transformations'])}; ordered: "
      f"{[o['parameters'] for o in parsed.get('parameter_dependencies', {}).get('ordered', [])]}")
print(f"  files: {run / 'fileManager.txt'}, {run / 'config.toml'}, {sset / lpi_name}")
