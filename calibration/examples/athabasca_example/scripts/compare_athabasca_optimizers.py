"""Compare Athabasca calibration runs (async DDS, SPOTPY DDS / SCE-UA, Pymoo GA) on the same J().

Usage (from my_work_athabasca/):
  python scripts/compare_athabasca_optimizers.py \
      "Async DDS=results/asyncdds_mpi_21217337" \
      "SPOTPY DDS=results/spotpy/dds_21343300" \
      "SPOTPY SCE-UA=results/spotpy/sceua_21339837" \
      [--out results/comparison]

Each run folder is read as written by its driver:
  async DDS     evaluations.csv        (n, kge, elapsed_s, ...)
  SPOTPY, GA    evaluations_by_time.csv (run, hours, kge, best_so_far_kge)
Every model run that was executed is counted, in the order it finished; failed runs (KGE <= -1e5)
count as evaluations but never as the best. Writes comparison.md, comparison.csv and
best_so_far_comparison.png to --out, and prints the table."""
import argparse, csv, json, sys
from pathlib import Path

PENALTY_CUT = -1e5
LEVELS = (0.80, 0.85, 0.8748, 0.88)       # 0.8748 = SYMFLUENCE's own reported best
AT_EVALS = (1000, 2500, 5000)
AT_HOURS = (12, 24)

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help='"Label=path/to/result_folder"')
ap.add_argument("--out", default="results/comparison")
args = ap.parse_args()


def load(folder):
    folder = Path(folder)
    if (folder / "evaluations_by_time.csv").exists():
        with open(folder / "evaluations_by_time.csv") as f:
            rows = [(float(r["hours"]), float(r["kge"])) for r in csv.DictReader(f)]
    elif (folder / "evaluations.csv").exists():
        with open(folder / "evaluations.csv") as f:
            rows = [(float(r["elapsed_s"]) / 3600, float(r["kge"])) for r in csv.DictReader(f)]
    else:
        sys.exit(f"{folder}: neither evaluations_by_time.csv nor evaluations.csv found")
    rows.sort()
    summ = json.loads((folder / "summary.json").read_text()) if (folder / "summary.json").exists() else {}
    return [h for h, _ in rows], [k for _, k in rows], summ


def best_curve(kge):
    best, out = float("-inf"), []
    for k in kge:
        if k > PENALTY_CUT and k > best:
            best = k
        out.append(best)
    return out


def fmt(v, nd=4):
    return "-" if v is None or v == float("-inf") else (f"{v:.{nd}f}" if isinstance(v, float) else f"{v:,}")


results = []
for spec in args.runs:
    label, _, path = spec.partition("=")
    if not path:
        sys.exit(f'give each run as "Label=folder", got: {spec}')
    hours, kge, summ = load(path)
    best = best_curve(kge)
    n, wall = len(kge), (summ.get("elapsed_s", hours[-1] * 3600) / 3600 if hours else 0.0)
    failed = sum(k <= PENALTY_CUT for k in kge)
    first = lambda lv: next((i for i, b in enumerate(best) if b >= lv), None)
    r = {"run": label, "folder": path, "evaluations": n, "wall_h": round(wall, 2),
         "evals_per_h": round(n / wall, 1) if wall else None, "failed": failed,
         "failed_pct": round(100 * failed / n, 1) if n else None, "best_kge": best[-1] if best else None}
    for lv in LEVELS:
        i = first(lv)
        r[f"evals_to_{lv}"] = None if i is None else i + 1
        r[f"hours_to_{lv}"] = None if i is None else round(hours[i], 2)
    for e in AT_EVALS:
        r[f"best_at_{e}_evals"] = best[min(e, n) - 1] if n >= e else None
    for h in AT_HOURS:
        idx = [i for i, t in enumerate(hours) if t <= h]
        r[f"best_at_{h}_h"] = best[idx[-1]] if idx and hours[-1] >= h else None
    r["_hours"], r["_best"] = hours, best
    results.append(r)

out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
lines = ["| | " + " | ".join(r["run"] for r in results) + " |", "|---" * (len(results) + 1) + "|"]
row = lambda name, key, nd=4: lines.append(f"| {name} | " + " | ".join(fmt(r[key], nd) for r in results) + " |")
row("Model runs executed", "evaluations")
row("Wall time (h)", "wall_h", 1)
row("Runs per hour", "evals_per_h", 0)
row("Failed runs (%)", "failed_pct", 1)
row("**Best KGE (2016-2019)**", "best_kge")
for e in AT_EVALS:
    row(f"Best KGE after {e:,} runs", f"best_at_{e}_evals")
for h in AT_HOURS:
    row(f"Best KGE after {h} h", f"best_at_{h}_h")
for lv in LEVELS:
    row(f"Runs to reach {lv}", f"evals_to_{lv}")
    row(f"Hours to reach {lv}", f"hours_to_{lv}", 1)
table = "\n".join(lines)
(out / "comparison.md").write_text(table + "\n")
keys = [k for k in results[0] if not k.startswith("_")]
with open(out / "comparison.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in results:
        w.writerow({k: r[k] for k in keys})
print(table)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    lo = 1.0
    for r in results:
        pts = [(i + 1, h, b) for i, (h, b) in enumerate(zip(r["_hours"], r["_best"])) if b > PENALTY_CUT]
        if not pts:
            continue
        lo = min(lo, pts[0][2])
        lab = f'{r["run"]} ({fmt(r["best_kge"])})'
        axes[0].step([p[0] for p in pts], [p[2] for p in pts], where="post", lw=2, label=lab)
        axes[1].step([p[1] for p in pts], [p[2] for p in pts], where="post", lw=2, label=lab)
    for ax, xl in zip(axes, ("Model runs", "Wall-clock hours")):
        ax.axhline(0.8748, ls=":", color="grey", lw=1)
        ax.set_xlabel(xl); ax.grid(alpha=0.3)
    axes[0].set_ylabel("Best KGE so far (2016-2019)")
    axes[0].set_ylim(max(0.0, lo - 0.02), min(1.0, max(r["best_kge"] or 0 for r in results) + 0.03))
    axes[1].legend(loc="lower right", title="Best KGE")
    axes[0].set_title("Same number of runs", fontsize=10); axes[1].set_title("Same wall time", fontsize=10)
    fig.suptitle("Athabasca (69 HRUs), J() develop build, 49 workers on 1 node; dotted: 0.8748 (SYMFLUENCE)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "best_so_far_comparison.png", dpi=150)
    print(f"\nFigure: {out / 'best_so_far_comparison.png'}")
except Exception as e:
    print(f"(figure skipped: {e})")
