"""Plot KGE per iteration for a SPOTPY run, next to SYMFLUENCE's iterations (works mid-run).
Usage: python plot_iterations.py [run_folder_name] [evaluations_per_iteration]
       e.g. python plot_iterations.py dds_full 50"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ross = Path(__file__).resolve().parent.parent
run = sys.argv[1] if len(sys.argv) > 1 else "dds_full"
batch = int(sys.argv[2]) if len(sys.argv) > 2 else 50
folder = ross / "results" / run
d = pd.read_csv(folder / "evaluations.csv")
d["best"] = d["kge"].cummax()
n = len(d)

# SPOTPY iterations: 0 = initial guess, then one per `batch` evaluations (last one may be partial)
ends = [1] + list(range(1 + batch, n + 1, batch))
if ends[-1] != n:
    ends.append(n)
spot = pd.DataFrame({"iteration": range(len(ends)), "evaluations": ends,
                     "best_kge": [d["best"].iloc[e - 1] for e in ends]})
sym = pd.read_csv(ross / "anvil-files" / "run_1_parallel_iteration_results.csv")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8))
ax1.plot(sym["iteration"], sym["score"], "o-", ms=3, label="SYMFLUENCE AsyncDDS (50 evaluations per iteration)")
ax1.plot(spot["iteration"], spot["best_kge"], "s-", ms=5,
         label=f"SPOTPY DDS ({batch} evaluations per iteration, {n} so far)")
ax1.set_xlabel("Iteration"); ax1.set_ylabel("Best KGE so far")
ax1.set_ylim(-0.2, 0.95); ax1.set_xlim(-0.5, 21); ax1.grid(alpha=0.3); ax1.legend(loc="lower right")
ax1.set_title("Calibration: best KGE per iteration (calibration period 2011–2015)")

ax2.scatter(d["n"], d["kge"].clip(lower=-1), s=6, alpha=0.4, label="each evaluation (values below −1 shown at −1)")
ax2.plot(d["n"], d["best"], color="C1", lw=2, label="best so far")
ax2.axhline(0.9053, color="gray", ls="--", lw=1, label="SYMFLUENCE final (0.9053)")
ax2.set_xlabel("SPOTPY evaluation number"); ax2.set_ylabel("KGE")
ax2.set_ylim(-1.05, 1); ax2.grid(alpha=0.3); ax2.legend(loc="lower right")
ax2.set_title(f"SPOTPY DDS: every evaluation ({(d.kge <= -1e5).sum()} failed runs)")

fig.tight_layout()
out = folder / "kge_per_iteration.png"
fig.savefig(out, dpi=150)
spot.to_csv(folder / "kge_per_iteration.csv", index=False)
print(f"Saved {out}")
print(spot.tail(5).to_string(index=False))
