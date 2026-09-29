"""Export a SPOTPY run in the same format as SYMFLUENCE's run_1_parallel_iteration_results.csv.
Usage: python export_iteration_results.py [run_folder_name] [evaluations_per_iteration]
       e.g. python export_iteration_results.py dds_full 50"""
import os, sys
from datetime import datetime, timedelta
from pathlib import Path
import pandas as pd

ross = Path(__file__).resolve().parent.parent
run = sys.argv[1] if len(sys.argv) > 1 else "dds_full"
batch = int(sys.argv[2]) if len(sys.argv) > 2 else 50
folder = ross / "results" / run
d = pd.read_csv(folder / "evaluations.csv")
params = [c for c in d.columns if c not in ("n", "kge", "seconds")]

# timestamps: run start (file creation time) plus cumulative run time
st = os.stat(folder / "evaluations.csv")
start = datetime.fromtimestamp(getattr(st, "st_birthtime", st.st_mtime - d["seconds"].sum()))
d["timestamp"] = [start + timedelta(seconds=s) for s in d["seconds"].cumsum()]
d["failed"] = d["kge"] <= -1e5

rows, ends = [], [1] + list(range(1 + batch, len(d) + 1, batch))
if ends[-1] != len(d):
    ends.append(len(d))                        # partial last iteration while the run is going
for it, end in enumerate(ends):
    part = d.iloc[:end]
    b = part["kge"].idxmax()
    batch_part = d.iloc[(ends[it - 1] if it else 0):end]
    row = {"iteration": it, "score": part.loc[b, "kge"],
           "timestamp": part["timestamp"].iloc[-1].isoformat()}
    row.update({p: f"[{part.loc[b, p]:.8g}]" for p in params})
    row["crash_count"] = int(batch_part["failed"].sum())
    row["crash_rate"] = round(batch_part["failed"].mean(), 4)
    row["evaluations"] = end
    rows.append(row)

out = folder / f"{run}_parallel_iteration_results.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"wrote {out}  ({len(rows)} rows, {len(d)} evaluations, {batch} per iteration)")
print(pd.DataFrame(rows)[["iteration", "evaluations", "score", "crash_count"]].tail(5).to_string(index=False))
