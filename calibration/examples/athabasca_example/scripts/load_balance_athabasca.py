#!/usr/bin/env python3
"""Load balancing of the Athabasca calibrations: parameter drain fraction and idle worker fraction.
  drain = (T_end - T_last_dispatch) / (T_end - T_first_dispatch)
  idle  = 1 - sum(run durations) / (n_workers * (T_end - T_first_dispatch))
Reads <folder>/evaluations.csv (async DDS / converted Fortran: elapsed_s; Pymoo GA: hours, generation).
"""
import argparse, pathlib, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(folder):
    f = pathlib.Path(folder) / "evaluations.csv"
    ranks = sorted(pathlib.Path(folder).glob("evaluations_rank_*.csv"))
    if f.exists():
        d = pd.read_csv(f)
    elif ranks:   # SPOTPY: one file per MPI rank, end_time in Unix seconds
        parts = [pd.read_csv(r).assign(worker=int(r.stem.rsplit("_", 1)[1])) for r in ranks]
        d = pd.concat([x for x in parts if len(x)], ignore_index=True)
        e = d["end_time"].astype(float)
        d["elapsed_s"] = e - (e - d["seconds"].astype(float)).min()
    else:
        print(f"  skipped: no evaluations.csv or evaluations_rank_*.csv in {folder}", file=sys.stderr); return None
    if "elapsed_s" in d: end = d["elapsed_s"].astype(float)
    elif "hours" in d: end = d["hours"].astype(float) * 3600.0
    else:
        print(f"  skipped: {f} has neither elapsed_s nor hours", file=sys.stderr); return None
    out = pd.DataFrame({"worker": d["worker"].astype(int),
                        "run_s": d["seconds"].astype(float).clip(lower=0.0), "end_s": end})
    out["start_s"] = (out["end_s"] - out["run_s"]).clip(lower=0.0)
    if "generation" in d: out["generation"] = d["generation"].astype(int)
    return out


def metrics(d, n_workers):
    t0, t1 = d.start_s.min(), d.end_s.max(); wall = t1 - t0
    m = {"Model runs": len(d), "Workers": n_workers,
         "Evaluations per worker": len(d) / n_workers,
         "Wall time from first dispatch (h)": wall / 3600,
         "Time before first dispatch (h)": t0 / 3600,
         "Median model run (min)": d.run_s[d.run_s > 60].median() / 60,
         "Last run (min)": (t1 - d.start_s.max()) / 60,
         "Drain fraction (%)": 100 * (t1 - d.start_s.max()) / wall,
         "Idle worker fraction (%)": 100 * (1 - d.run_s.sum() / (n_workers * wall))}
    gens = None
    if "generation" in d:
        gens = d.groupby("generation").agg(runs=("run_s", "size"), first_dispatch_s=("start_s", "min"),
                                           end_s=("end_s", "max"), busy_s=("run_s", "sum"),
                                           median_run_s=("run_s", "median"), slowest_run_s=("run_s", "max"))
        gens["wall_s"] = gens.end_s - gens.first_dispatch_s
        gens["idle_pct"] = 100 * (1 - gens.busy_s / (n_workers * gens.wall_s))
        gens["gap_before_s"] = gens.first_dispatch_s - gens.end_s.shift(1)
        m["Idle within generations (%)"] = 100 * (1 - gens.busy_s.sum() / (n_workers * gens.wall_s.sum()))
    return m, gens


def busy_curve(d):
    t = np.r_[d.start_s.values, d.end_s.values]
    step = np.r_[np.ones(len(d)), -np.ones(len(d))]
    o = np.lexsort((step, t))
    return t[o], np.cumsum(step[o])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help='"Label=folder" with evaluations.csv')
    ap.add_argument("--workers", type=int, default=49)
    ap.add_argument("--out", default="results/load_balance")
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows, data = {}, {}
    for spec in a.runs:
        label, folder = spec.split("=", 1)
        print(f"{label}: {folder}", file=sys.stderr)
        d = load(folder)
        if d is None: continue
        if d.worker.nunique() != a.workers:
            print(f"  note: {d.worker.nunique()} distinct workers in the file, using --workers {a.workers}", file=sys.stderr)
        m, gens = metrics(d, a.workers)
        rows[label], data[label] = m, d
        if gens is not None:
            name = "generations_" + label.lower().replace(" ", "_") + ".csv"
            gens.round(2).to_csv(out / name)
            print(f"  per-generation idle: {out / name}", file=sys.stderr)
    if not rows: sys.exit("no runs could be read")
    t = pd.DataFrame(rows)
    fmt = {"Model runs": "{:,.0f}", "Workers": "{:.0f}"}
    shown = t.apply(lambda r: r.map(lambda v: "-" if pd.isna(v) else fmt.get(r.name, "{:.1f}").format(v)), axis=1)
    md = "| | " + " | ".join(shown.columns) + " |\n|---|" + "---|" * len(shown.columns) + "\n"
    md += "".join(f"| {i} | " + " | ".join(shown.loc[i]) + " |\n" for i in shown.index)
    print(md); (out / "load_balance.md").write_text(md); t.to_csv(out / "load_balance.csv")
    labels = list(rows)
    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.45, wspace=0.25)
    for k, (key, title, sub) in enumerate([
            ("Drain fraction (%)", "Parameter drain fraction", "Time after final parameter dispatch"),
            ("Idle worker fraction (%)", "Total idle worker fraction", "Idle worker capacity over the full calibration")]):
        ax = fig.add_subplot(gs[0, k]); v = [rows[l][key] for l in labels]
        bars = ax.bar(labels, v, color="0.35", width=0.6)
        ax.bar_label(bars, [f"{x:.1f}%" for x in v], padding=3, fontsize=9)
        ax.set_title(f"{title}\n", loc="left", fontsize=12)
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=9)
        ax.set_ylabel("Fraction of runtime (%)" if k == 0 else "Fraction of worker capacity (%)")
        ax.set_ylim(0, max(v) * 1.25 + 0.5); ax.grid(axis="y", color="0.9"); ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=9)
    ax = fig.add_subplot(gs[1, :])
    for l in labels:
        tt, b = busy_curve(data[l])
        ax.step(tt / 3600, 100 * b / a.workers, where="post", lw=0.7,
                label=f"{l} (idle {rows[l]['Idle worker fraction (%)']:.1f}%)")
    ax.set_xlabel("Hours since start of calibration"); ax.set_ylabel("Workers running a model (%)")
    ax.set_ylim(0, 105); ax.set_title("Busy workers over time", loc="left", fontsize=12)
    ax.grid(color="0.9"); ax.legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Athabasca calibrations: load balancing ({a.workers} workers)", x=0.06, ha="left", fontsize=14, weight="bold")
    fig.savefig(out / "load_balance.png", dpi=150, bbox_inches="tight")
    print(f"Figure: {out / 'load_balance.png'}", file=sys.stderr)


if __name__ == "__main__":
    main()
