"""SPOTPY calibration of Athabasca (69 HRUs) with J(params): NCAR/summa develop build, noTS config,
period-end observations. Same model setup as calibrate_athabasca_asyncdds_mpi_develop.py; only the
optimizer differs (SPOTPY instead of the async DDS driver).

Usage (from summa-spotpy/):
  mpirun -np N python -u my_work_athabasca/scripts/spotpy/calibrate_athabasca_spotpy.py <sceua|dds|lhs> <repetitions> [ngs] [max_loop_inc] [tag]
    ngs           SCE-UA complexes (default N-1 = one per worker; in SPOTPY each worker evolves one complex)
    max_loop_inc  SCE-UA only: 0 = no limit (default); 1 = burn-in only (smoke test)
    dds           SPOTPY DDS (Tolson & Shoemaker 2007), r = 0.2, starting from the default parameters;
                  with parallel="mpi" each free worker gets a candidate perturbed from the current best

Rank 0 is SPOTPY's MPI master (for DDS it also runs the starting, default-parameter evaluation); ranks 1..N-1 run one J() evaluation at a time, each in
its own child process (scripts/eval_one_develop.py), so a fatal SUMMA error ends only that evaluation.
Writes only to results/spotpy/<algo>_<tag>/ and output/spotpy/<algo>_<tag>/rank_NNN/."""
import csv, json, math, os, re, subprocess, sys, time
from pathlib import Path
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()

basin = Path(__file__).resolve().parents[2]                 # .../my_work_athabasca
ALGO = sys.argv[1] if len(sys.argv) > 1 else "sceua"
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
NGS = int(sys.argv[3]) if len(sys.argv) > 3 and int(sys.argv[3]) > 0 else max(size - 1, 2)
MAXLOOP = int(sys.argv[4]) if len(sys.argv) > 4 else 0
TAG = sys.argv[5] if len(sys.argv) > 5 else os.environ.get("SLURM_JOB_ID", "local")
assert ALGO in ("sceua", "dds", "lhs"), "algorithm must be sceua, dds or lhs"
MAXIMIZE = ALGO in ("lhs", "dds")                           # SCE-UA minimises -> objective = -KGE
EVAL_TIMEOUT, PENALTY, SEED = 2 * 3600, -1e6, 42            # as in the async DDS driver

out = basin / "results" / "spotpy" / f"{ALGO}_{TAG}"
work = basin / "output" / "spotpy" / f"{ALGO}_{TAG}" / f"rank_{rank:03d}"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

with open(basin / "scripts" / "athabasca_params.csv") as f:
    rows = list(csv.DictReader(f))
NAMES = [r["name"] for r in rows]

# ---- per-rank setup, identical to the async DDS develop driver. Rank 0 gets it too: SPOTPY's DDS
#      evaluates the starting (default) parameter set on the master itself.
if True:
    fm = [f"outputPath '{work}/'" if l.strip().startswith("outputPath") else ("outputControlFile 'outputControl_develop.txt'" if l.strip().startswith("outputControlFile") else l)
          for l in (basin / "settings" / "SUMMA" / "fileManager.txt").read_text().splitlines()]
    (work / "fileManager.txt").write_text("\n".join(fm) + "\n")
    os.environ["SUMMA_MASTER"] = str(work / "fileManager.txt")
    topo = work / "topology.nc"
    if not topo.exists():
        os.symlink(basin / "settings" / "mizuRoute" / "topology.nc", topo)
    cfg = (basin / "settings" / "athabasca_config_develop_noTS.toml").read_text()
    cfg = re.sub(r'^hfabric_path\s*=.*$', f'hfabric_path = "{work}/"', cfg, flags=re.M)
    (work / "athabasca_config.toml").write_text(cfg)
    os.environ["SUMMA_CONFIG"] = str(work / "athabasca_config.toml")
    os.chdir(work)
    console = open(work / "summa_console.log", "a", buffering=1)   # SUMMA's daily printing goes here
    if rank > 0:
        os.dup2(console.fileno(), 1); os.dup2(console.fileno(), 2)

import numpy as np
import spotpy
# SPOTPY reseeds numpy from random_state when the sampler is created, so the seed is passed there (below)
evaluator = basin / "scripts" / "eval_one_develop.py"
evlog = open(out / f"evaluations_rank_{rank:03d}.csv", "w", buffering=1)
evlog.write("n,kge,seconds,end_time," + ",".join(NAMES) + "\n")

class AthabascaSetup:
    def __init__(self):
        self.params = [spotpy.parameter.Uniform(r["name"], low=float(r["low"]), high=float(r["high"]),
                                                optguess=float(r["default"])) for r in rows]
        self.n = 0
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    def simulation(self, x):
        x, t0 = [float(v) for v in x], time.time()
        self.n += 1
        p, kge = dict(zip(NAMES, x)), PENALTY
        if p["heightCanopyTop"] > p["heightCanopyBottom"]:
            res = work / f"kge_{self.n}.txt"
            try:
                cp = subprocess.run([sys.executable, "-u", str(evaluator), str(res), json.dumps(NAMES), json.dumps(x)],
                                    cwd=work, timeout=EVAL_TIMEOUT, stdout=console, stderr=subprocess.STDOUT)
                if cp.returncode == 0 and res.exists():
                    kge = float(res.read_text())
                else:
                    print(f"[driver] rank {rank} evaluation {self.n} crashed (exit code {cp.returncode}); scored as failed", file=console, flush=True)
            except subprocess.TimeoutExpired:
                print(f"[driver] rank {rank} evaluation {self.n} exceeded {EVAL_TIMEOUT} s; scored as failed", file=console, flush=True)
            except Exception as e:
                print(f"[driver] rank {rank} evaluation {self.n} failed: {e}", file=console, flush=True)
            res.unlink(missing_ok=True)
        if not math.isfinite(kge):
            kge = PENALTY
        if evlog:
            evlog.write(f"{self.n},{kge},{time.time() - t0:.1f},{time.time():.0f}," + ",".join(map(str, x)) + "\n")
        return [kge]
    def evaluation(self):
        return [0.0]
    def objectivefunction(self, simulation, evaluation, params=None):
        return simulation[0] if MAXIMIZE else -simulation[0]

if rank == 0:
    print(f"SPOTPY {spotpy.__version__} {ALGO}: {REPS} repetitions, {size - 1} workers"
          + (f", ngs={NGS}, burn-in {NGS * (2 * len(NAMES) + 1)} runs" if ALGO == "sceua" else "")
          + (", r=0.2, starting from the default parameters" if ALGO == "dds" else "")
          + (", burn-in only" if (MAXLOOP == 1 and ALGO == "sceua") else "") + f", seed {SEED}", flush=True)
t_start = time.time()
kw = dict(dbname=str(out / f"{ALGO}_results"), dbformat="csv", parallel="mpi", save_sim=False,
          random_state=SEED + rank)   # rank 0 (master, proposes all candidates) uses SEED
if ALGO == "sceua":
    sampler = spotpy.algorithms.sceua(AthabascaSetup(), **kw)
    sampler.sample(REPS, ngs=NGS, max_loop_inc=(MAXLOOP or None))
elif ALGO == "dds":
    sampler = spotpy.algorithms.dds(AthabascaSetup(), r=0.2, **kw)
    sampler.sample(REPS, trials=1, x_initial=[float(r["default"]) for r in rows])
else:
    sampler = spotpy.algorithms.lhs(AthabascaSetup(), **kw)
    sampler.sample(REPS)

# ---- master only (SPOTPY's workers exit inside the sampler) ----
if rank == 0:
    elapsed = time.time() - t_start
    data = sampler.getdata()
    like = np.array(data["like1"], dtype=float)
    kge = like if MAXIMIZE else -like
    ok = kge > PENALTY / 10
    i = int(np.nanargmax(kge))
    # all runs actually executed (SPOTPY finishes each evolution loop and discards results past the
    # budget), in the order they finished, from the per-worker logs
    import glob
    ev = []
    for fn in sorted(glob.glob(str(out / "evaluations_rank_*.csv"))):
        with open(fn) as f:
            for r in csv.DictReader(f):
                ev.append((float(r["end_time"]), float(r["kge"]), float(r["seconds"])))
    ev.sort()
    t_end = np.array([e[0] for e in ev]); k_run = np.array([e[1] for e in ev])
    hours = (t_end - t_start) / 3600
    best_t = np.maximum.accumulate(np.where(k_run > PENALTY / 10, k_run, -np.inf))
    with open(out / "evaluations_by_time.csv", "w") as f:
        f.write("run,hours,kge,best_so_far_kge\n")
        for n, (h, k, b) in enumerate(zip(hours, k_run, best_t), 1):
            f.write(f"{n},{h:.3f},{k},{b}\n")
    def first_hours(level):
        idx = np.nonzero(best_t >= level)[0]
        return round(float(hours[idx[0]]), 2) if idx.size else None
    summary = {"algorithm": ALGO, "spotpy_version": spotpy.__version__, "processes": size, "workers": size - 1,
               "ngs": NGS if ALGO == "sceua" else None, "repetitions_requested": REPS,
               "evaluations_saved": int(len(kge)), "runs_executed": int(len(ev)),
               "failed_saved": int((~ok).sum()), "elapsed_s": round(elapsed, 1),
               "runs_per_hour": round(len(ev) / elapsed * 3600, 1), "seed": SEED,
               "best_kge": float(kge[i]), "hours_to_kge": {str(l): first_hours(l) for l in (0.80, 0.85, 0.8748, 0.88)},
               "best_params": {n: float(data["par" + n][i]) for n in NAMES}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print("SUMMARY " + json.dumps({k: v for k, v in summary.items() if k != "best_params"}), flush=True)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
        fin = np.isfinite(best_t)
        ax.plot(hours[fin], best_t[fin], lw=2, label=f"SPOTPY {ALGO.upper()} ({len(ev)} runs)")
        ax.axhline(0.8806, ls="--", color="grey", label="Async DDS develop best (0.8806, 28.4 h)")
        ax.set_xlabel("Wall-clock hours"); ax.set_ylabel("Best KGE so far (2016–2019)")
        ax.set_ylim(max(0.0, float(best_t[fin].min()) - 0.05), 1.0)
        ax.grid(alpha=0.3); ax.legend(loc="lower right"); fig.tight_layout()
        fig.savefig(out / "best_so_far.png", dpi=150)
    except Exception as e:
        print(f"(plot skipped: {e})", flush=True)
