#!/usr/bin/env python3
"""Parameter-level load balancing, in the style of the "Balancing parameter trials" slide.
Each calibration is cut into consecutive blocks of 1,000 runs (dispatch order, measured durations),
replayed with asynchronous scheduling on 192/95/47/23 workers = 5.2/10.5/21.3/43.5 evaluations per worker.
Last column: the actual runs as executed.
"""
import argparse, heapq, pathlib, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKERS = [192, 95, 47, 23]
COLORS = ["#2a6fdb", "#e07b00", "#1b9e77", "#c2185b", "#7b5cc4", "#8c6d31", "#555555"]


def load(folder):
    folder = pathlib.Path(folder)
    f = folder / "evaluations.csv"
    ranks = sorted(folder.glob("evaluations_rank_*.csv"))
    if f.exists():
        d = pd.read_csv(f)
    elif ranks:
        parts = [pd.read_csv(r).assign(worker=int(r.stem.rsplit("_", 1)[1])) for r in ranks]
        d = pd.concat([x for x in parts if len(x)], ignore_index=True)
        for c in ("seconds", "end_time"): d[c] = pd.to_numeric(d[c])
        d["elapsed_s"] = d.end_time - (d.end_time - d.seconds).min()
    else:
        print(f"  skipped: no evaluations in {folder}", file=sys.stderr); return None
    end = d["elapsed_s"].astype(float) if "elapsed_s" in d else d["hours"].astype(float) * 3600.0
    run = pd.to_numeric(d["seconds"]).astype(float).clip(lower=0.0)
    out = pd.DataFrame({"worker": d["worker"].astype(int), "run_s": run, "end_s": end})
    out["start_s"] = (out.end_s - out.run_s).clip(lower=0.0)
    return out.sort_values("start_s").reset_index(drop=True)


def actual(d, n_workers):
    t0, t1 = d.start_s.min(), d.end_s.max(); wall = t1 - t0
    return 100 * (t1 - d.start_s.max()) / wall, 100 * (1 - d.run_s.sum() / (n_workers * wall))


def replay(durations, n_workers):
    free = [0.0] * n_workers
    last_dispatch = end = 0.0
    for s in durations:
        t = heapq.heappop(free)
        last_dispatch = t
        heapq.heappush(free, t + s)
        end = max(end, t + s)
    return 100 * (end - last_dispatch) / end, 100 * (1 - sum(durations) / (n_workers * end))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help='"Label=folder"')
    ap.add_argument("--block", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=49)
    ap.add_argument("--out", default="results/load_balance")
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in a.runs:
        label, folder = spec.split("=", 1)
        print(f"{label}: {folder}", file=sys.stderr)
        d = load(folder)
        if d is None: continue
        dur = d.run_s.values
        for b in range(len(dur) // a.block):
            block = dur[b * a.block:(b + 1) * a.block]
            for w in WORKERS:
                dr, idle = replay(block, w)
                rows.append(dict(method=label, block=b + 1, workers=w, evals_per_worker=a.block / w,
                                 x=f"{a.block / w:.1f}", drain=dr, idle=idle, kind="simulated"))
        dr, idle = actual(d, a.workers)
        rows.append(dict(method=label, block=0, workers=a.workers, evals_per_worker=len(d) / a.workers,
                         x="actual", drain=dr, idle=idle, kind="actual"))
    if not rows: sys.exit("no runs could be read")
    t = pd.DataFrame(rows)
    t.round(3).to_csv(out / "load_balance_scaling.csv", index=False)
    sim = t[t.kind == "simulated"]
    summ = sim.groupby("workers").agg(evals_per_worker=("evals_per_worker", "first"), blocks=("drain", "size"),
                                      drain_median=("drain", "median"), idle_median=("idle", "median"),
                                      idle_max=("idle", "max")).sort_index(ascending=False)
    print(summ.round(2).to_string())
    print("\nActual runs (49 workers):")
    print(t[t.kind == "actual"][["method", "evals_per_worker", "drain", "idle"]].round(2).to_string(index=False))
    methods = list(dict.fromkeys(t.method))
    color = {m: COLORS[i % len(COLORS)] for i, m in enumerate(methods)}
    cats = [f"{a.block / w:.1f}" for w in WORKERS]
    rng = np.random.default_rng(0)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "0.15"})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    for ax, key, title, sub, ylab in [
            (axes[0], "drain", "Parameter drain fraction", "Time after final parameter dispatch", "Fraction of calibration runtime"),
            (axes[1], "idle", "Total idle worker fraction", "Idle worker capacity over the full calibration", "Fraction of worker capacity")]:
        for i, c in enumerate(cats):
            v = sim[sim.x == c]
            ax.boxplot(v[key], positions=[i], widths=0.32, showfliers=False, whis=(0, 100),
                       medianprops=dict(color="black", lw=1.4), boxprops=dict(lw=1.1), whiskerprops=dict(lw=1.1),
                       capprops=dict(lw=0))
            ax.scatter(i + rng.uniform(-0.07, 0.07, len(v)), v[key], s=14, alpha=0.45,
                       c=[color[m] for m in v.method], edgecolors="none", zorder=3)
            ax.scatter([i], [v[key].mean()], s=60, c="black", zorder=4)
        i = len(cats) + 0.35
        ax.axvline(len(cats) - 0.35, color="0.6", lw=0.8, ls=":")
        act = t[t.kind == "actual"]
        for m, val in zip(act.method, act[key]):
            ax.scatter(i + rng.uniform(-0.06, 0.06), val, s=46, c=color[m], edgecolors="black", lw=0.5, zorder=4)
        ax.set_xticks(list(range(len(cats))) + [i])
        ax.set_xticklabels(cats + [f"Actual runs\n({a.workers} workers, ~100+)"])
        ax.set_xlim(-0.6, i + 0.6)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))
        ax.grid(color="0.9", lw=0.8); ax.set_axisbelow(True)
        ax.set_xlabel("Parameter evaluations per worker", fontsize=12)
        ax.set_ylabel(ylab, fontsize=12)
        ax.set_title(title, loc="left", fontsize=15, pad=24)
        ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=11.5)
    handles = [plt.Line2D([], [], marker="o", ls="", color=color[m], label=m) for m in methods]
    fig.legend(handles=handles, loc="lower center", ncol=len(methods), frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    fig.text(0.01, 0.995, "Parameter-level load balancing: Athabasca calibrations", fontsize=15, va="top")
    fig.text(0.01, 0.95, f"Each point is one block of {a.block:,} model runs with their measured run times, replayed "
             "with asynchronous scheduling; last column: the calibrations as actually run", fontsize=11, va="top")
    fig.tight_layout(rect=(0, 0.06, 1, 0.9))
    for ext in ("png", "pdf"):
        fig.savefig(out / f"load_balance_scaling.{ext}", dpi=170, bbox_inches="tight")
    print(f"\nFigure: {out / 'load_balance_scaling.png'} (and .pdf)", file=sys.stderr)


if __name__ == "__main__":
    main()
