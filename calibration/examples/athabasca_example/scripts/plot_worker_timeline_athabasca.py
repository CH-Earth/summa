#!/usr/bin/env python3
"""Worker activity through time, coloured by KGE: one panel per calibration run.
Each bar is one model run on one worker, from dispatch to completion; white gaps are idle time.
Options: --sharey (same time axis for all panels), --ncols 3, --kge-min -1
"""
import argparse, pathlib, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

FAILED = -1e5   # penalty scores (-1e6) and anything below are failed runs


def load(folder):
    folder = pathlib.Path(folder)
    f = folder / "evaluations.csv"
    ranks = sorted(folder.glob("evaluations_rank_*.csv"))
    if f.exists():
        d = pd.read_csv(f)
    elif ranks:   # SPOTPY: one file per MPI rank, end_time in Unix seconds
        parts = [pd.read_csv(r).assign(worker=int(r.stem.rsplit("_", 1)[1])) for r in ranks]
        d = pd.concat([x for x in parts if len(x)], ignore_index=True)
        e = d["end_time"].astype(float)
        d["elapsed_s"] = e - (e - d["seconds"].astype(float)).min()
    else:
        print(f"  skipped: no evaluations in {folder}", file=sys.stderr); return None
    if "elapsed_s" in d: end = d["elapsed_s"].astype(float)
    elif "hours" in d: end = d["hours"].astype(float) * 3600.0
    else:
        print(f"  skipped: {folder} has neither elapsed_s nor hours", file=sys.stderr); return None
    run = d["seconds"].astype(float).clip(lower=0.0)
    out = pd.DataFrame({"worker": d["worker"].astype(int), "run_s": run, "end_s": end,
                        "kge": d["kge"].astype(float)})
    out["start_s"] = (out.end_s - out.run_s).clip(lower=0.0)
    out[["start_s", "end_s"]] -= out.start_s.min()
    return out


def panel(ax, d, cmap, norm, gap_s):
    good = d.kge > FAILED
    w = 0.42
    for mask, colors in ((good, cmap(norm(d.kge[good].clip(lower=norm.vmin)))), (~good, "#d9d9d9")):
        x, s, e = d.worker[mask].values, d.start_s[mask].values / 3600, (d.end_s[mask].values - gap_s) / 3600
        e = np.maximum(e, s)
        verts = [((xi - w, si), (xi + w, si), (xi + w, ei), (xi - w, ei)) for xi, si, ei in zip(x, s, e)]
        ax.add_collection(PolyCollection(verts, facecolors=colors, edgecolors="none"))
    ax.set_xlim(d.worker.min() - 0.7, d.worker.max() + 0.7)
    ax.set_ylim(0, d.end_s.max() / 3600 * 1.01)
    ax.set_facecolor("white")
    ax.grid(color="0.92", lw=0.5); ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help='"Label=folder"')
    ap.add_argument("--out", default="results/load_balance")
    ap.add_argument("--ncols", type=int, default=3)
    ap.add_argument("--sharey", action="store_true")
    ap.add_argument("--kge-min", type=float, default=-1.0)
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    runs = []
    for spec in a.runs:
        label, folder = spec.split("=", 1)
        print(f"{label}: {folder}", file=sys.stderr)
        d = load(folder)
        if d is not None: runs.append((label, d))
    if not runs: sys.exit("no runs could be read")
    ncols = min(a.ncols, len(runs)); nrows = -(-len(runs) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows), squeeze=False, sharey=a.sharey)
    cmap = plt.get_cmap("viridis"); norm = matplotlib.colors.Normalize(vmin=a.kge_min, vmax=1.0)
    for ax, (label, d) in zip(axes.flat, runs):
        panel(ax, d, cmap, norm, gap_s=0.06 * d.run_s[d.run_s > 60].median())
        n_fail = int((d.kge <= FAILED).sum())
        ax.set_title(f"{label}\n{len(d):,} runs, {d.end_s.max() / 3600:.1f} h, best KGE {d.kge.max():.3f}"
                     + (f", {n_fail} failed (grey)" if n_fail else ""), fontsize=10, loc="left")
    for ax in axes.flat[len(runs):]: ax.set_visible(False)
    if a.sharey:
        top = max(d.end_s.max() for _, d in runs) / 3600 * 1.01
        for ax in axes.flat: ax.set_ylim(0, top)
    for ax in axes[:, 0]: ax.set_ylabel("Hours since first dispatch")
    for ax in axes.flat[:len(runs)]: ax.set_xlabel("Worker rank")
    fig.suptitle("Athabasca calibrations: worker activity through time, coloured by KGE",
                 x=0.01, ha="left", fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 0.93, 0.97))
    cax = fig.add_axes((0.945, 0.3, 0.012, 0.4))
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(f"KGE (values below {a.kge_min:g} shown at the bottom colour)")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"worker_timeline.{ext}", dpi=170)
    print(f"Figure: {out / 'worker_timeline.png'} (and .pdf)", file=sys.stderr)


if __name__ == "__main__":
    main()
