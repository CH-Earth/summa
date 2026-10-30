"""Pymoo genetic-algorithm calibration of Athabasca (69 HRUs) with J(params): NCAR/summa develop build,
noTS config, period-end observations. Same model setup as calibrate_athabasca_asyncdds_mpi_develop.py and
scripts/spotpy/calibrate_athabasca_spotpy.py; only the optimizer differs (Pymoo GA).

Usage (from summa-spotpy/):
  mpirun -np N python -u my_work_athabasca/scripts/pymoo/calibrate_athabasca_pymoo.py \
      [--evals 5000] [--pop 98] [--offspring N-1] [--seed 42] [--timeout 7200] [--tag TAG] [--resume]

  GA         : pymoo GA (real-coded; SBX crossover + polynomial mutation, tournament selection,
               elitist survival), all parameters in their linear ranges from athabasca_params.csv,
               like the async DDS and SPOTPY runs.
  start      : the default parameter set + (pop - 1) Latin hypercube samples (DDS starts from the defaults).
  generation : 'offspring' new parameter sets (default N-1 = one per worker, so one generation is one
               wave across all workers). Budget: pop + offspring * generations <= evals
               (pop 98, offspring 49, evals 5000 -> 1 + 100 generations = 4,998 evaluations).
  checkpoint : results/pymoo/ga_<tag>/checkpoint.pkl after every generation; rerun with the same --tag
               and --resume to continue a run that hit its walltime (the interrupted generation is redone).

Rank 0 runs the GA and hands out parameter sets to whichever worker is free; it runs no SUMMA itself.
Ranks 1..N-1 run one J() evaluation at a time, each in its own child process (scripts/eval_one_develop.py)
with a timeout, so a fatal SUMMA error or a hung run ends only that evaluation (KGE = -1e6).
Writes only to results/pymoo/ga_<tag>/ and output/pymoo/ga_<tag>/rank_NNN/."""
import argparse, csv, json, math, os, re, subprocess, sys, time
from pathlib import Path
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()
TAG_TASK, TAG_RESULT, TAG_STOP = 1, 2, 3

ap = argparse.ArgumentParser()
ap.add_argument("--evals", type=int, default=5000, help="evaluation budget (default 5000)")
ap.add_argument("--pop", type=int, default=98, help="population size (default 98 = two waves of 49)")
ap.add_argument("--offspring", type=int, default=0, help="new sets per generation (default: one per worker)")
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--timeout", type=int, default=2 * 3600, help="seconds per evaluation (default 7200, as async DDS)")
ap.add_argument("--tag", default=os.environ.get("SLURM_JOB_ID", "local"))
ap.add_argument("--resume", action="store_true", help="continue from results/pymoo/ga_<tag>/checkpoint.pkl")
args = ap.parse_args()
NWORK = size - 1
assert NWORK >= 1, "needs at least 2 MPI processes (1 coordinator + workers)"
NOFF = args.offspring if args.offspring > 0 else NWORK
assert args.pop >= 2 and args.evals >= args.pop, "need evals >= pop >= 2"
NGEN = 1 + (args.evals - args.pop) // NOFF                 # generation 1 = the initial population
PENALTY = -1e6                                            # as in the async DDS and SPOTPY drivers

basin = Path(__file__).resolve().parents[2]               # .../my_work_athabasca
out = basin / "results" / "pymoo" / f"ga_{args.tag}"
work = basin / "output" / "pymoo" / f"ga_{args.tag}" / f"rank_{rank:03d}"
with open(basin / "scripts" / "athabasca_params.csv") as f:
    rows = list(csv.DictReader(f))
NAMES = [r["name"] for r in rows]


# =====================================================================================================
# Workers: same per-rank SUMMA setup and J() call as the async DDS develop driver and the SPOTPY driver
# =====================================================================================================
def worker():
    work.mkdir(parents=True, exist_ok=True)
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
    os.dup2(console.fileno(), 1); os.dup2(console.fileno(), 2)
    evaluator = basin / "scripts" / "eval_one_develop.py"
    st = MPI.Status()
    while True:
        msg = comm.recv(source=0, tag=MPI.ANY_TAG, status=st)
        if st.Get_tag() == TAG_STOP:
            break
        eid, x = msg
        t0, kge, status = time.time(), PENALTY, "ok"
        p = dict(zip(NAMES, x))
        if not p["heightCanopyTop"] > p["heightCanopyBottom"]:
            status = "canopy_top<=bottom"
        else:
            res = work / f"kge_{eid}.txt"
            try:
                cp = subprocess.run([sys.executable, "-u", str(evaluator), str(res), json.dumps(NAMES), json.dumps(x)],
                                    cwd=work, timeout=args.timeout, stdout=console, stderr=subprocess.STDOUT)
                if cp.returncode == 0 and res.exists():
                    kge = float(res.read_text())
                else:
                    status = f"crashed({cp.returncode})"
            except subprocess.TimeoutExpired:
                status = "timeout"
            except Exception as e:
                status = f"error({type(e).__name__})"
            res.unlink(missing_ok=True)
            if status == "ok" and not math.isfinite(kge):
                kge, status = PENALTY, "nonfinite"
            if status == "ok" and kge <= PENALTY / 10:
                status = "J_returned_penalty"     # J() itself reported a failed run
            if status != "ok":
                print(f"[driver] rank {rank} evaluation {eid}: {status}; scored as failed", file=console, flush=True)
        comm.send((eid, kge, status, time.time() - t0, time.time()), dest=0, tag=TAG_RESULT)


# =====================================================================================================
# Coordinator (rank 0): the GA itself
# =====================================================================================================
def coordinator():
    import dill
    import numpy as np
    import pymoo
    from pymoo.algorithms.soo.nonconvex.ga import GA
    from pymoo.core.evaluator import Evaluator
    from pymoo.core.problem import Problem
    from pymoo.operators.sampling.lhs import LHS
    from pymoo.problems.static import StaticProblem

    out.mkdir(parents=True, exist_ok=True)
    xl = np.array([float(r["low"]) for r in rows]); xu = np.array([float(r["high"]) for r in rows])
    x_default = np.array([float(r["default"]) for r in rows])
    problem = Problem(n_var=len(NAMES), n_obj=1, xl=xl, xu=xu)
    ckpt, evfile, genfile = out / "checkpoint.pkl", out / "evaluations.csv", out / "generations.csv"

    if args.resume:
        state = dill.loads(ckpt.read_bytes())
        algorithm, elapsed_before, n_eval = state["algorithm"], state["elapsed_s"], state["n_eval"]
        gen_done = algorithm.n_gen - 1                     # pymoo counts the next generation after tell()
        # drop rows from the generation that was interrupted; it is evaluated again
        for fn in (evfile, genfile):
            with open(fn) as f:
                lines = f.readlines()
            keep = [lines[0]] + [l for l in lines[1:] if int(l.split(",")[1 if fn == evfile else 0]) <= gen_done]
            fn.write_text("".join(keep))
        print(f"Resuming {out.name} after generation {gen_done} ({n_eval} evaluations, "
              f"{elapsed_before / 3600:.2f} h so far)", flush=True)
    else:
        assert not ckpt.exists(), f"{ckpt} exists: use --resume, or a new --tag"
        rng = np.random.default_rng(args.seed)
        X0 = LHS().do(problem, args.pop - 1, random_state=rng).get("X")
        X0 = np.vstack([x_default, X0])                     # the defaults are individual 1
        algorithm = GA(pop_size=args.pop, n_offsprings=NOFF, sampling=X0, eliminate_duplicates=True)
        algorithm.setup(problem, termination=("n_gen", NGEN), seed=args.seed)
        elapsed_before, n_eval = 0.0, 0
        evfile.write_text("n,generation,worker,kge,status,seconds,hours," + ",".join(NAMES) + "\n")
        genfile.write_text("generation,evaluations,hours,generation_best_kge,best_so_far_kge,failed,"
                           "slowest_run_s,median_run_s\n")
        for name, src in (("athabasca_params.csv", basin / "scripts" / "athabasca_params.csv"),
                          ("athabasca_config_develop_noTS.toml", basin / "settings" / "athabasca_config_develop_noTS.toml")):
            (out / name).write_text(src.read_text())        # record the settings used

    print(f"Pymoo {pymoo.__version__} GA: pop {args.pop}, {NOFF} offspring per generation, {NGEN} generations "
          f"= {args.pop + (NGEN - 1) * NOFF} evaluations (budget {args.evals}), {NWORK} workers, "
          f"timeout {args.timeout} s, seed {args.seed}; start = defaults + {args.pop - 1} LHS", flush=True)
    t_start = time.time()
    hours_now = lambda: (elapsed_before + time.time() - t_start) / 3600
    best_so_far = [float(-algorithm.opt.get("F").min())] if algorithm.opt is not None else [-np.inf]
    st = MPI.Status()

    def evaluate(X, gen):
        """Evaluate all rows of X: hand each to whichever worker is free. Returns KGE per row."""
        nonlocal n_eval
        todo = list(range(len(X)))[::-1]
        kge, secs = np.full(len(X), np.nan), np.zeros(len(X))
        busy = 0
        with open(evfile, "a", buffering=1) as log:
            for w in range(1, NWORK + 1):
                if not todo:
                    break
                i = todo.pop(); comm.send((i, [float(v) for v in X[i]]), dest=w, tag=TAG_TASK); busy += 1
            while busy:
                i, k, status, s, _ = comm.recv(source=MPI.ANY_SOURCE, tag=TAG_RESULT, status=st)
                w = st.Get_source(); busy -= 1
                kge[i], secs[i] = k, s
                n_eval += 1
                log.write(f"{n_eval},{gen},{w},{k},{status},{s:.1f},{hours_now():.4f},"
                          + ",".join(repr(float(v)) for v in X[i]) + "\n")
                if todo:
                    j = todo.pop(); comm.send((j, [float(v) for v in X[j]]), dest=w, tag=TAG_TASK); busy += 1
        return kge, secs

    while algorithm.has_next():
        pop = algorithm.ask()
        gen = algorithm.n_gen if algorithm.n_gen is not None else 1
        X = pop.get("X")
        kge, secs = evaluate(X, gen)
        Evaluator().eval(StaticProblem(problem, F=-kge.reshape(-1, 1)), pop)   # pymoo minimises -KGE
        algorithm.tell(infills=pop)
        ok = kge > PENALTY / 10
        best_so_far.append(max(best_so_far[-1], float(np.max(kge))))
        with open(genfile, "a") as f:
            f.write(f"{gen},{n_eval},{hours_now():.4f},{np.max(kge)},{best_so_far[-1]},{int((~ok).sum())},"
                    f"{secs.max():.0f},{np.median(secs):.0f}\n")
        print(f"GA generation {gen}/{NGEN}: {n_eval} evaluations, {hours_now():.2f} h, "
              f"best KGE {best_so_far[-1]:.4f} (this generation {np.max(kge):.4f}), failed {int((~ok).sum())}, "
              f"slowest run {secs.max() / 60:.1f} min", flush=True)
        tmp = ckpt.with_suffix(".tmp")
        tmp.write_bytes(dill.dumps({"algorithm": algorithm, "n_eval": n_eval,
                                    "elapsed_s": elapsed_before + time.time() - t_start}))
        os.replace(tmp, ckpt)                                # atomic: a kill never leaves half a checkpoint

    for w in range(1, NWORK + 1):
        comm.send(None, dest=w, tag=TAG_STOP)

    # ---- results, in the same form as the SPOTPY driver ----
    elapsed = elapsed_before + time.time() - t_start
    with open(evfile) as f:
        ev = list(csv.DictReader(f))
    hours = np.array([float(r["hours"]) for r in ev]); k_run = np.array([float(r["kge"]) for r in ev])
    best_t = np.maximum.accumulate(np.where(k_run > PENALTY / 10, k_run, -np.inf))
    with open(out / "evaluations_by_time.csv", "w") as f:
        f.write("run,hours,kge,best_so_far_kge\n")
        for n, (h, k, b) in enumerate(zip(hours, k_run, best_t), 1):
            f.write(f"{n},{h:.3f},{k},{b}\n")
    i = int(np.argmax(np.where(k_run > PENALTY / 10, k_run, -np.inf)))

    def first(level, axis):
        idx = np.nonzero(best_t >= level)[0]
        return None if not idx.size else (round(float(hours[idx[0]]), 2) if axis == "h" else int(idx[0]) + 1)

    summary = {"algorithm": "pymoo GA", "pymoo_version": pymoo.__version__, "processes": size, "workers": NWORK,
               "population": args.pop, "offspring_per_generation": NOFF, "generations_incl_initial": int(algorithm.n_gen - 1),
               "evaluation_budget": args.evals, "evaluations": len(ev),
               "failed": int((k_run <= PENALTY / 10).sum()),
               "failed_by_reason": {s: sum(r["status"] == s for r in ev) for s in sorted({r["status"] for r in ev}) if s != "ok"},
               "elapsed_s": round(elapsed, 1), "runs_per_hour": round(len(ev) / elapsed * 3600, 1),
               "eval_timeout_s": args.timeout, "seed": args.seed, "start": "defaults + LHS",
               "default_params_kge": next((float(r["kge"]) for r in ev
                                           if all(abs(float(r[n]) - d) <= 1e-12 * max(1.0, abs(d)) for n, d in zip(NAMES, x_default))), None),
               "best_kge": float(k_run[i]), "best_found_at_evaluation": i + 1,
               "hours_to_kge": {str(l): first(l, "h") for l in (0.80, 0.85, 0.8748, 0.88)},
               "evaluations_to_kge": {str(l): first(l, "n") for l in (0.80, 0.85, 0.8748, 0.88)},
               "best_params": {n: float(ev[i][n]) for n in NAMES}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print("SUMMARY " + json.dumps({k: v for k, v in summary.items() if k != "best_params"}), flush=True)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        fin = np.isfinite(best_t)
        for ax, xv, xlab in ((axes[0], np.arange(1, len(ev) + 1), "Evaluations"), (axes[1], hours, "Wall-clock hours")):
            ax.plot(xv[fin], best_t[fin], lw=2, label=f"Pymoo GA ({len(ev)} runs)")
            ax.axhline(0.8806, ls="--", color="grey", label="Async DDS develop best (0.8806, 28.4 h)")
            ax.set_xlabel(xlab); ax.grid(alpha=0.3)
        axes[0].set_ylabel("Best KGE so far (2016–2019)")
        axes[0].set_ylim(max(0.0, float(best_t[fin].min()) - 0.05), 1.0)
        axes[1].legend(loc="lower right"); fig.tight_layout()
        fig.savefig(out / "best_so_far.png", dpi=150)
    except Exception as e:
        print(f"(plot skipped: {e})", flush=True)


if rank == 0:
    try:
        coordinator()
    except BaseException:
        import traceback; traceback.print_exc(); sys.stdout.flush(); sys.stderr.flush()
        comm.Abort(1)       # otherwise the workers would wait for work until the walltime
else:
    worker()
