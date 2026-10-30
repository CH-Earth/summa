"""[develop build, noTS config, period-end obs] Asynchronous parallel DDS calibration of Athabasca, v2.

Same as calibrate_athabasca_asyncdds_mpi_develop.py (v1: same J(), workers, child-process isolation,
2 h timeout, -1e6 penalty, seed 42, r = 0.2), plus the search features of the SUMMA develop Fortran DDS
(build/source/calibration/summa_parameter_sampling.f90 and parameter_search.f90):

  - log-space search for parameters marked "log" in the settings file ([transformations]);
  - ordered parameter constraints with a minimum gap ([[ordered]]), checked BEFORE a model run:
    candidates that break them are redrawn and cost no evaluation;
  - (option) out-of-bound steps rejected and redrawn, as in the Fortran DDS, instead of reflected;
  - (option) random first wave, as in the Fortran DDS, instead of starting from the defaults.
  - the default parameter set runs in the first wave together with 48 DDS candidates (v1 ran it alone).

Usage (from summa-spotpy/):
  mpirun -np N python -u my_work_athabasca/scripts/calibrate_athabasca_asyncdds_mpi_develop_v2.py \
      [max_evals=5000] [r=0.2] [tag] [--settings FILE] [--bounds reflect|reject] [--start defaults|random]
  defaults: settings = scripts/athabasca_dds_v2_settings.toml, bounds = reflect, start = defaults
  Check the settings first, on a login node:  python .../calibrate_athabasca_asyncdds_mpi_develop_v2.py --check
  Writes results/asyncdds_mpi_v2_<tag>/ (same files as v1, plus search_settings.json) and
  output/asyncdds_v2_<tag>/rank_NNN/."""
import argparse, csv, json, math, os, shutil, sys, time
from pathlib import Path
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()
TAG_WORK, TAG_RESULT, TAG_STOP = 1, 2, 0

basin = Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser()
ap.add_argument("max_evals", nargs="?", type=int, default=5000)
ap.add_argument("r", nargs="?", type=float, default=0.2)
ap.add_argument("tag", nargs="?", default=os.environ.get("SLURM_JOB_ID", "local"))
ap.add_argument("--settings", default=str(basin / "scripts" / "athabasca_dds_v2_settings.toml"))
ap.add_argument("--bounds", choices=("reflect", "reject"), default="reflect",
                help="reflect: Tolson & Shoemaker (2007), as v1; reject: redraw, as the Fortran DDS")
ap.add_argument("--start", choices=("defaults", "random"), default="defaults",
                help="defaults: as v1; random: random first wave, as the Fortran DDS")
ap.add_argument("--check", action="store_true",
                help="print the resolved search settings and exit (no MPI job needed)")
args = ap.parse_args()
M, R, TAG = args.max_evals, args.r, args.tag
EVAL_TIMEOUT = 2 * 3600   # seconds per evaluation
BATCH, PENALTY, SEED, MAXTRY = 50, -1e6, 42, 10000

out = basin / "results" / f"asyncdds_mpi_v2_{TAG}"
work = basin / "output" / f"asyncdds_v2_{TAG}" / f"rank_{rank:03d}"
if not args.check:
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

with open(basin / "scripts" / "athabasca_params.csv") as f:
    rows = list(csv.DictReader(f))
NAMES = [r["name"] for r in rows]
LOW = [float(r["low"]) for r in rows]
HIGH = [float(r["high"]) for r in rows]
X0 = [float(r["default"]) for r in rows]

# =============================== worker (unchanged from v1) ===============================
if rank > 0:
    fm = [f"outputPath '{work}/'" if l.strip().startswith("outputPath") else ("outputControlFile 'outputControl_develop.txt'" if l.strip().startswith("outputControlFile") else l)
          for l in (basin / "settings" / "SUMMA" / "fileManager.txt").read_text().splitlines()]
    (work / "fileManager.txt").write_text("\n".join(fm) + "\n")
    os.environ["SUMMA_MASTER"] = str(work / "fileManager.txt")
    import re
    topo = work / "topology.nc"
    if not topo.exists():
        os.symlink(basin / "settings" / "mizuRoute" / "topology.nc", topo)
    cfg = (basin / "settings" / "athabasca_config_develop_noTS.toml").read_text()
    cfg = re.sub(r'^hfabric_path\s*=.*$', f'hfabric_path = "{work}/"', cfg, flags=re.M)
    (work / "athabasca_config.toml").write_text(cfg)
    os.environ["SUMMA_CONFIG"] = str(work / "athabasca_config.toml")
    os.chdir(work)
    fd = os.open(str(work / "summa_console.log"), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.dup2(fd, 1); os.dup2(fd, 2)                                   # SUMMA's daily printing goes to a per-rank log
    # Each evaluation runs J() in its own child process, so a fatal SUMMA/mizuRoute error or a hung
    # run ends only that evaluation (scored PENALTY) and this worker carries on.
    import subprocess
    evaluator = basin / "scripts" / "eval_one_develop.py"
    status = MPI.Status()
    while True:
        msg = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        if status.Get_tag() == TAG_STOP:
            break
        n, x = msg
        t0 = time.time()
        p = dict(zip(NAMES, x))
        kge = PENALTY
        if p["heightCanopyTop"] > p["heightCanopyBottom"]:                # safety net; the master already checks
            res = work / f"kge_{n}.txt"
            try:
                cp = subprocess.run([sys.executable, "-u", str(evaluator), str(res), json.dumps(NAMES), json.dumps(list(x))],
                                    cwd=work, timeout=EVAL_TIMEOUT)
                if cp.returncode == 0 and res.exists():
                    kge = float(res.read_text())
                else:
                    print(f"[driver] evaluation {n} crashed (exit code {cp.returncode}); scored as failed", flush=True)
            except subprocess.TimeoutExpired:
                print(f"[driver] evaluation {n} exceeded {EVAL_TIMEOUT} s; scored as failed", flush=True)
            except Exception as e:
                print(f"[driver] evaluation {n} failed: {e}", flush=True)
            res.unlink(missing_ok=True)
        if not math.isfinite(kge):
            kge = PENALTY
        comm.send((n, kge, time.time() - t0), dest=0, tag=TAG_RESULT)
    sys.exit(0)

# =============================== master (rank 0) ===============================
import numpy as np


def stop_all(msg):
    print("ERROR: " + msg, flush=True)
    sys.exit(1) if size == 1 else comm.Abort(1)


try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
settings_file = Path(args.settings)
settings = tomllib.loads(settings_file.read_text())

# ---- transformations (natural log, as in parameter_search.f90) ----
TRANS = settings.get("transformations", {})
for name, t in TRANS.items():
    if name not in NAMES:
        stop_all(f"[transformations] {name} is not one of the calibrated parameters")
    if t not in ("log", "none"):
        stop_all(f"[transformations] {name}: unsupported transformation {t!r}")
    if t == "log" and LOW[NAMES.index(name)] <= 0:
        stop_all(f"[transformations] {name}: log needs a positive lower bound")
is_log = np.array([TRANS.get(n, "none") == "log" for n in NAMES])
low, high = np.array(LOW), np.array(HIGH)
s_low = np.where(is_log, np.log(np.where(is_log, low, 1.0)), low)       # search-space bounds
s_high = np.where(is_log, np.log(np.where(is_log, high, 1.0)), high)
s_span = s_high - s_low
to_search = lambda x: np.where(is_log, np.log(np.where(is_log, np.maximum(x, 1e-300), 1.0)), x)
to_phys = lambda s: np.where(is_log, np.exp(np.where(is_log, s, 0.0)), s)

# ---- fixed (non-calibrated) parameters used in ordered constraints ----
def read_localparaminfo():
    vals = {}
    f = basin / "settings" / "SUMMA" / "localParamInfo.txt"
    if f.exists():
        for line in f.read_text().splitlines():
            parts = [p.strip() for p in line.split("|")]
            if len(parts) == 4 and parts[0] and not parts[0].startswith(("!", "'", "*")):
                try:
                    vals[parts[0]] = {"value": float(parts[1].replace("d", "e")),
                                      "lower": float(parts[2].replace("d", "e")),
                                      "upper": float(parts[3].replace("d", "e")), "source": "localParamInfo.txt"}
                except ValueError:
                    pass
    return vals

lpi = read_localparaminfo()
fixed = {}
for name, v in settings.get("fixed", {}).items():
    fixed[name] = {"value": float(v["value"]), "lower": float(v.get("lower", v["value"])),
                   "upper": float(v.get("upper", v["value"])), "source": "settings [fixed]"}

ORDERED = []   # list of (names, gap, [("search", index) | ("fixed", value)])
for k, oc in enumerate(settings.get("ordered", []), 1):
    names, gf = list(oc["parameters"]), float(oc.get("gap_fraction", 0.0))
    if len(names) < 2 or gf < 0:
        stop_all(f"[[ordered]] #{k}: needs >= 2 parameters and gap_fraction >= 0")
    members, lows, highs = [], [], []
    for n in names:
        if n in NAMES:
            j = NAMES.index(n); members.append(("search", j)); lows.append(LOW[j]); highs.append(HIGH[j])
        else:
            info = fixed.get(n) or lpi.get(n)
            if info is None:
                stop_all(f"[[ordered]] #{k}: {n} is not calibrated and not found in [fixed] or localParamInfo.txt")
            fixed[n] = info
            members.append(("fixed", info["value"])); lows.append(info["lower"]); highs.append(info["upper"])
    gap = gf * (highs[-1] - lows[0])           # as initialize_constraints in parameter_search.f90
    ORDERED.append((names, gap, members))
    # feasibility check, as in parameter_search.f90
    prev = None
    for n, (kind, v), lo, hi in zip(names, members, lows, highs):
        need = lo if prev is None else max(lo, prev + gap)
        if kind == "search":
            if need > hi:
                stop_all(f"[[ordered]] #{k} is infeasible at {n}: needs >= {need:.4g}, upper bound {hi:.4g}")
            prev = need
        else:
            if v < need:
                stop_all(f"[[ordered]] #{k}: fixed {n} = {v:.4g} is below the required {need:.4g}")
            prev = v


def satisfies_ordered(x):
    for names, gap, members in ORDERED:
        prev = None
        for kind, v in members:
            val = x[v] if kind == "search" else v
            if prev is not None and val < prev + gap:
                return False
            prev = val
    return True


x0 = np.array(X0)
if not satisfies_ordered(x0):
    stop_all("the default parameter set breaks an [[ordered]] constraint; adjust gap_fraction or the defaults")

rng = np.random.default_rng(SEED)
redraws = {"bounds": 0, "ordered": 0}


def dds_candidate(xbest, i):
    """DDS perturbation of the current best in search space; i = candidate number (1-based)."""
    p = 1.0 - math.log(i) / math.log(M) if i > 1 else 1.0
    sb = to_search(xbest)
    for _ in range(MAXTRY):
        sel = rng.random(len(sb)) < p
        if not sel.any():
            sel[rng.integers(len(sb))] = True
        s = sb.copy()
        rejected = False
        for d in np.where(sel)[0]:
            v = s[d] + R * s_span[d] * rng.standard_normal()
            if args.bounds == "reject":
                if v < s_low[d] or v > s_high[d]:
                    rejected = True
                    break
            elif v < s_low[d]:                         # reflection (Tolson & Shoemaker 2007), as v1
                v = s_low[d] + (s_low[d] - v)
                if v > s_high[d]:
                    v = s_low[d]
            elif v > s_high[d]:
                v = s_high[d] - (v - s_high[d])
                if v < s_low[d]:
                    v = s_high[d]
            s[d] = v
        if rejected:
            redraws["bounds"] += 1
            continue
        x = np.clip(to_phys(s), low, high)             # clip only guards exp/log round-off
        if satisfies_ordered(x):
            return x
        redraws["ordered"] += 1
    stop_all(f"could not generate a DDS candidate satisfying the constraints in {MAXTRY} tries")


def random_candidate():
    for _ in range(MAXTRY):
        x = np.clip(to_phys(s_low + rng.random(len(s_low)) * s_span), low, high)
        if satisfies_ordered(x):
            return x
        redraws["ordered"] += 1
    stop_all(f"could not generate a random set satisfying the constraints in {MAXTRY} tries")


if args.check:
    print(f"Settings file: {settings_file}")
    print(f"log-space parameters: {[n for n, l in zip(NAMES, is_log) if l] or 'none'}")
    for names, gap, members in ORDERED:
        shown = [f"{n}={v:.4g} (fixed, from {fixed[n]['source']})" if k == "fixed" else n
                 for n, (k, v) in zip(names, members)]
        print(f"ordered: {' < '.join(shown)}, minimum gap {gap:.4g}")
    print(f"default parameter set satisfies all constraints: True")
    print(f"bounds rule: {args.bounds}, start: {args.start}")
    sys.exit(0)
assert size >= 2, "need at least 2 MPI processes (1 master + 1 worker)"

search_settings = {
    "version": "v2", "settings_file": str(settings_file), "bounds": args.bounds, "start": args.start,
    "log_parameters": [n for n, l in zip(NAMES, is_log) if l],
    "ordered": [{"parameters": names, "min_gap": gap} for names, gap, _ in ORDERED],
    "fixed_parameters": fixed,
}
shutil.copy(basin / "settings" / "athabasca_config_develop_noTS.toml", out)
shutil.copy(basin / "scripts" / "athabasca_params.csv", out)
shutil.copy(settings_file, out)
(out / "search_settings.json").write_text(json.dumps(search_settings, indent=2))
evf = open(out / "evaluations.csv", "w", buffering=1)
evf.write("n,worker,kge,seconds,elapsed_s,best_so_far," + ",".join(NAMES) + "\n")
itf = open(out / "iterations_best.csv", "w", buffering=1)
itf.write("iteration,evaluations,elapsed_s,batch_best_kge,failed_runs,best_so_far_kge," + ",".join(NAMES) + "\n")
print(f"Async DDS v2: {M} evaluations, r={R}, {size - 1} workers, seed {SEED}, bounds={args.bounds}, "
      f"start={args.start}", flush=True)
print(f"  log-space parameters: {search_settings['log_parameters'] or 'none'}", flush=True)
for names, gap, members in ORDERED:
    shown = [f"{n}={v:.4g} (fixed)" if k == "fixed" else n for n, (k, v) in zip(names, members)]
    print(f"  ordered: {' < '.join(shown)}, minimum gap {gap:.4g}", flush=True)

t_start, status = time.time(), MPI.Status()
pending, state = {}, {"n_sent": 0, "n_done": 0, "batch_best": -np.inf, "batch_fail": 0}
xbest, fbest = x0.copy(), -np.inf


def send(w, x):
    state["n_sent"] += 1
    pending[state["n_sent"]] = x
    comm.send((state["n_sent"], [float(v) for v in x]), dest=w, tag=TAG_WORK)


def receive_and_record():
    global xbest, fbest
    (n, kge, sec) = comm.recv(source=MPI.ANY_SOURCE, tag=TAG_RESULT, status=status)
    w = status.Get_source()
    x = pending.pop(n)
    state["n_done"] += 1
    if kge >= fbest:                                  # DDS keeps the new point if not worse
        xbest, fbest = x, kge
    state["batch_best"] = max(state["batch_best"], kge)
    state["batch_fail"] += int(kge <= PENALTY / 10)
    el = time.time() - t_start
    evf.write(f"{n},{w},{kge},{sec:.1f},{el:.1f},{fbest}," + ",".join(repr(float(v)) for v in x) + "\n")
    nd = state["n_done"]
    if nd % BATCH == 0 or nd == M:
        it = math.ceil(nd / BATCH)
        itf.write(f"{it},{nd},{el:.1f},{state['batch_best']},{state['batch_fail']},{fbest}," +
                  ",".join(repr(float(v)) for v in xbest) + "\n")
        print(f"AsyncDDS-v2 {it}/{math.ceil(M / BATCH)} | evals {nd} | best KGE {fbest:.4f} | "
              f"batch best {state['batch_best']:.4f} | failed {state['batch_fail']} | "
              f"redraws {redraws['bounds']}+{redraws['ordered']} | elapsed {el / 60:.1f} min", flush=True)
        state["batch_best"], state["batch_fail"] = -np.inf, 0
    return w


# 1) first wave: all workers start at once.
#    start=defaults: the default set + DDS candidates around it (v1 ran the defaults alone first)
#    start=random:   random sets, as the Fortran DDS does
for w in range(1, size):
    if state["n_sent"] >= M:
        break
    if args.start == "defaults":
        send(w, x0.copy() if w == 1 else dds_candidate(x0, state["n_sent"] + 1))
    else:
        send(w, random_candidate())
# 2) refill each worker as soon as it returns, with a candidate around the current best
while state["n_done"] < M:
    w = receive_and_record()
    if state["n_sent"] < M:
        send(w, dds_candidate(xbest, state["n_sent"] + 1))
for w in range(1, size):
    comm.send(None, dest=w, tag=TAG_STOP)

elapsed = time.time() - t_start
summary = {"algorithm": "async DDS v2", "processes": size, "workers": size - 1, "evaluations": M,
           "r": R, "seed": SEED, "bounds": args.bounds, "start": args.start,
           "log_parameters": search_settings["log_parameters"],
           "redraws_out_of_bounds": redraws["bounds"], "redraws_ordered": redraws["ordered"],
           "elapsed_s": round(elapsed, 1), "evals_per_hour": round(M / elapsed * 3600, 1), "best_kge": fbest,
           "best_params": dict(zip(NAMES, map(float, xbest)))}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print("SUMMARY " + json.dumps({k: v for k, v in summary.items() if k != "best_params"}), flush=True)

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = pd.read_csv(out / "iterations_best.csv")
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(df["iteration"], df["best_so_far_kge"], lw=2, label="Best KGE so far (async DDS v2, J())")
    ax.scatter(df["iteration"], df["batch_best_kge"].clip(lower=-1), s=12, alpha=0.5, label="Best KGE in iteration")
    ax.axhline(0.8806, ls="--", color="grey", label="Async DDS v1 develop best (0.8806)")
    ax.set_xlabel(f"Iteration ({BATCH} evaluations each)")
    ax.set_ylabel("KGE (calibration period, see athabasca_config.toml)")
    ax.set_title(f"Athabasca async DDS v2: {M} evaluations, {size - 1} workers, {elapsed / 3600:.2f} h")
    ax.set_ylim(max(-1, float(df["batch_best_kge"].min()) - 0.05), 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "kge_per_iteration.png", dpi=150)
    print(f"Plot: {out / 'kge_per_iteration.png'}", flush=True)
except Exception as e:
    print(f"Plot skipped: {e}", flush=True)
