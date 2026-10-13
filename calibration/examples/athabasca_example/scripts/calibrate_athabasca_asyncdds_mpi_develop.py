"""[develop build, noTS config, period-end obs] Asynchronous parallel DDS calibration of Athabasca (69 HRUs) with J(params) (MPI master-worker).
Usage: mpirun -np N python calibrate_athabasca_asyncdds_mpi.py [max_evals=5000] [r=0.2] [tag]
Rank 0 runs DDS (Tolson & Shoemaker 2007) and hands out candidates; ranks 1..N-1 run SUMMA.
Each time a worker returns, the best solution is updated and that worker immediately receives a
new candidate perturbed from the current best, so no worker waits for a batch to finish.
One 'iteration' = 50 completed evaluations (as in SYMFLUENCE's batches of 50)."""
import csv, json, math, os, shutil, sys, time
from pathlib import Path
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()
TAG_WORK, TAG_RESULT, TAG_STOP = 1, 2, 0

basin = Path(__file__).resolve().parent.parent
M = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
R = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2
TAG = sys.argv[3] if len(sys.argv) > 3 else os.environ.get("SLURM_JOB_ID", "local")
EVAL_TIMEOUT = 2 * 3600   # seconds per evaluation
BATCH, PENALTY, SEED, SYMFLUENCE_KGE = 50, -1e6, 42, None   # set to SYMFLUENCE's Athabasca best KGE for the plot

out = basin / "results" / f"asyncdds_mpi_{TAG}"
work = basin / "output" / f"asyncdds_{TAG}" / f"rank_{rank:03d}"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

with open(basin / "scripts" / "athabasca_params.csv") as f:
    rows = list(csv.DictReader(f))
NAMES = [r["name"] for r in rows]
LOW = [float(r["low"]) for r in rows]
HIGH = [float(r["high"]) for r in rows]
X0 = [float(r["default"]) for r in rows]

# =============================== worker ===============================
if rank > 0:
    fm = [f"outputPath '{work}/'" if l.strip().startswith("outputPath") else ("outputControlFile 'outputControl_develop.txt'" if l.strip().startswith("outputControlFile") else l)
          for l in (basin / "settings" / "SUMMA" / "fileManager.txt").read_text().splitlines()]
    (work / "fileManager.txt").write_text("\n".join(fm) + "\n")
    os.environ["SUMMA_MASTER"] = str(work / "fileManager.txt")
    # per-rank config: mizuRoute reads the river network from, and writes its subset file into, the rank folder
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
    # Each evaluation runs J() in its own child process (scripts/eval_one.py). A fatal Fortran
    # error (error stop) in SUMMA/mizuRoute then ends only that child: the evaluation is scored
    # as failed (PENALTY), as SYMFLUENCE does, and this worker carries on. Runs longer than
    # EVAL_TIMEOUT are stopped and scored as failed too (SYMFLUENCE: SUMMA_TIMEOUT 7200 s).
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
        if p["heightCanopyTop"] > p["heightCanopyBottom"]:
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
assert size >= 2, "need at least 2 MPI processes (1 master + 1 worker)"
rng = np.random.default_rng(SEED)
low, high = np.array(LOW), np.array(HIGH)
span = high - low

def dds_candidate(xbest, i):
    """DDS perturbation of the current best; i = candidate number (1-based)."""
    p = 1.0 - math.log(i) / math.log(M) if i > 1 else 1.0
    sel = rng.random(len(xbest)) < p
    if not sel.any():
        sel[rng.integers(len(xbest))] = True
    x = xbest.copy()
    for d in np.where(sel)[0]:
        v = x[d] + R * span[d] * rng.standard_normal()
        if v < low[d]:
            v = low[d] + (low[d] - v)
            if v > high[d]:
                v = low[d]
        elif v > high[d]:
            v = high[d] - (v - high[d])
            if v < low[d]:
                v = high[d]
        x[d] = v
    return x

shutil.copy(basin / "settings" / "athabasca_config_develop_noTS.toml", out)
shutil.copy(basin / "scripts" / "athabasca_params.csv", out)
evf = open(out / "evaluations.csv", "w", buffering=1)
evf.write("n,worker,kge,seconds,elapsed_s,best_so_far," + ",".join(NAMES) + "\n")
itf = open(out / "iterations_best.csv", "w", buffering=1)
itf.write("iteration,evaluations,elapsed_s,batch_best_kge,failed_runs,best_so_far_kge," + ",".join(NAMES) + "\n")
print(f"Async DDS: {M} evaluations, r={R}, {size - 1} workers, seed {SEED}", flush=True)

t_start, status = time.time(), MPI.Status()
pending, state = {}, {"n_sent": 0, "n_done": 0, "batch_best": -np.inf, "batch_fail": 0}
xbest, fbest = np.array(X0), -np.inf

def send(w, x):
    state["n_sent"] += 1
    pending[state["n_sent"]] = x
    comm.send((state["n_sent"], x.tolist()), dest=w, tag=TAG_WORK)

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
    evf.write(f"{n},{w},{kge},{sec:.1f},{el:.1f},{fbest}," + ",".join(map(str, x)) + "\n")
    nd = state["n_done"]
    if nd % BATCH == 0 or nd == M:
        it = math.ceil(nd / BATCH)
        itf.write(f"{it},{nd},{el:.1f},{state['batch_best']},{state['batch_fail']},{fbest}," +
                  ",".join(map(str, xbest)) + "\n")
        print(f"AsyncDDS {it}/{math.ceil(M / BATCH)} | evals {nd} | best KGE {fbest:.4f} | "
              f"batch best {state['batch_best']:.4f} | failed {state['batch_fail']} | "
              f"elapsed {el / 60:.1f} min", flush=True)
        state["batch_best"], state["batch_fail"] = -np.inf, 0
    return w

# 1) initial guess (default parameters), as SYMFLUENCE injects it
send(1, np.array(X0))
receive_and_record()
print(f"Initial guess KGE = {fbest:.4f}", flush=True)

# 2) fill every worker, then refill each one as soon as it returns
for w in range(1, size):
    if state["n_sent"] < M:
        send(w, dds_candidate(xbest, state["n_sent"] + 1))
while state["n_done"] < M:
    w = receive_and_record()
    if state["n_sent"] < M:
        send(w, dds_candidate(xbest, state["n_sent"] + 1))
for w in range(1, size):
    comm.send(None, dest=w, tag=TAG_STOP)

elapsed = time.time() - t_start
summary = {"algorithm": "async DDS", "processes": size, "workers": size - 1, "evaluations": M,
           "r": R, "seed": SEED, "elapsed_s": round(elapsed, 1),
           "evals_per_hour": round(M / elapsed * 3600, 1), "best_kge": fbest,
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
    ax.plot(df["iteration"], df["best_so_far_kge"], lw=2, label="Best KGE so far (async DDS, J())")
    ax.scatter(df["iteration"], df["batch_best_kge"].clip(lower=-1), s=12, alpha=0.5, label="Best KGE in iteration")
    if SYMFLUENCE_KGE is not None:
        ax.axhline(SYMFLUENCE_KGE, ls="--", color="grey", label=f"SYMFLUENCE async DDS best ({SYMFLUENCE_KGE})")
    ax.set_xlabel(f"Iteration ({BATCH} evaluations each)")
    ax.set_ylabel("KGE (calibration period, see athabasca_config.toml)")
    ax.set_title(f"Athabasca async DDS on Anvil: {M} evaluations, {size - 1} workers, {elapsed / 3600:.2f} h")
    ax.set_ylim(max(-1, float(df["batch_best_kge"].min()) - 0.05), 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "kge_per_iteration.png", dpi=150)
    print(f"Plot: {out / 'kge_per_iteration.png'}", flush=True)
except Exception as e:
    print(f"Plot skipped: {e}", flush=True)
