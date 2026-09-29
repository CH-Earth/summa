"""Progress report for a SPOTPY calibration run.
Usage: python check_progress.py [run_folder_name] [total_evaluations]
       e.g. python check_progress.py dds_full 1000"""
import subprocess, sys
from datetime import datetime, timedelta
from pathlib import Path
import pandas as pd

ross = Path(__file__).resolve().parent.parent
run = sys.argv[1] if len(sys.argv) > 1 else "dds_full"
total = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
folder = ross / "results" / run
log = folder / "evaluations.csv"

running = "calibrate_ross_spotpy.py" in subprocess.run(["ps", "aux"], capture_output=True, text=True).stdout
print(f"Run: {run}   ({'RUNNING' if running else 'not running'})")
if not log.exists():
    sys.exit("No evaluations yet (evaluations.csv not found).")

d = pd.read_csv(log)
n = len(d)
if n == 0:
    sys.exit("Started, waiting for the first evaluation to finish.")
s = d["seconds"].tail(50).mean()
left = max(total - n, 0) * s
ib = d["kge"].idxmax()
print(f"Progress: {n}/{total} evaluations ({100*n/total:.0f}%), {s:.0f} s each (last 50)")
if running and left > 0:
    print(f"Time left: about {left/3600:.1f} h, finishing around {(datetime.now()+timedelta(seconds=left)):%a %H:%M}")
print(f"Best KGE: {d.kge.iloc[ib]:.4f} (evaluation {ib+1})   failed runs: {(d.kge <= -1e5).sum()}")

sym = ross / "anvil-files" / "run_1_parallel_iteration_results.csv"
if sym.exists():
    h = pd.read_csv(sym)
    ref = h.loc[h["iteration"] * 50 + 1 <= n, "score"].max()
    print(f"SYMFLUENCE best after {n} evaluations: {ref:.4f}   (final 0.9053)")

print("\nLast 5 evaluations:")
print(d[["n", "kge", "seconds"]].tail(5).to_string(index=False))

out = folder / "run.out"
if out.exists():
    errs = [l for l in out.read_text(errors="ignore").splitlines() if "Traceback" in l or "Error" in l]
    if errs:
        print("\nErrors in run.out:\n  " + "\n  ".join(errs[-3:]))
