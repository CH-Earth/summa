import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

nodes = [1, 2, 4, 8]
wall = [12.333, 6.343, 3.670, 2.257]          # walltime_h from scaling_tests/nodes_*/run_info.txt
ideal = [wall[0] / n for n in nodes]
eff = [100 * wall[0] / (n * t) for n, t in zip(nodes, wall)]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
fig.subplots_adjust(left=0.07, right=0.98, top=0.72, bottom=0.12, wspace=0.28)
fig.text(0.035, 0.965, "SUMMA model-instance strong scaling", fontsize=13)
fig.text(0.035, 0.925, "Athabasca (69 HRUs) single-basin calibration; 5000 SUMMA evaluations; "
         "async DDS, 127 workers per node", fontsize=10.5)

def style(ax):
    for s in ax.spines.values():
        s.set_color("black"); s.set_linewidth(0.8)
    ax.grid(True, color="#E5E5E5", linewidth=0.8); ax.set_axisbelow(True)
    ax.tick_params(colors="#4D4D4D", length=3)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(FixedLocator(nodes))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0.9, 8.9)
    ax.set_xlabel("Compute nodes", fontsize=11, labelpad=8)

# --- left: wall time
ax = axes[0]
ax.plot(nodes, ideal, "--", color="black", lw=1.8, dashes=(5, 3), label="Ideal")
ax.plot(nodes, wall, "-", color="black", lw=1.8, label="Observed")
ax.plot(nodes, wall, "o", color="black", ms=8)
ax.set_yscale("log", base=2)
ax.yaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_ylim(1.3, 14.5)
ax.set_ylabel("Wall-clock time (hours)", fontsize=11)
style(ax)
for n, t in zip(nodes, wall):
    off = (-44, 8) if n == 8 else (8, 6)
    ax.annotate(f"{t:.2f} h", (n, t), textcoords="offset points", xytext=off, fontsize=9, color="#333333")
ax.set_title("Strong scaling\n", loc="left", fontsize=12, pad=26)
ax.text(0, 1.115, "Athabasca; 5000 evaluations; 127 workers per node", transform=ax.transAxes, fontsize=10)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=9.5, handlelength=2.5)

# --- right: efficiency
ax = axes[1]
ax.axhline(100, ls="--", color="black", lw=1.8, dashes=(5, 3))
ax.plot(nodes, eff, "-", color="black", lw=1.8)
ax.plot(nodes, eff, "o", color="black", ms=8)
ax.set_ylim(60, 106)
ax.yaxis.set_major_locator(FixedLocator([60, 70, 80, 90, 100]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
ax.set_ylabel("Parallel efficiency", fontsize=11)
style(ax)
for n, e in zip(nodes, eff):
    if n > 1:
        ax.annotate(f"{e:.0f}%", (n, e), textcoords="offset points", xytext=(-26, -14), fontsize=9, color="#333333")
ax.set_title("Strong-scaling efficiency\n", loc="left", fontsize=12, pad=26)
ax.text(0, 1.115, "Efficiency relative to the 1-node run", transform=ax.transAxes, fontsize=10)

fig.savefig("athabasca_strong_scaling.pdf")
fig.savefig("athabasca_strong_scaling.png", dpi=200)
print([round(e, 1) for e in eff], [round(wall[0] / t, 2) for t in wall])
