"""Trace les sorties de advection_demo : instantanés + évolution des diagnostics.

Usage : python tools/plot_advection.py <dossier_des_sorties> [figure.png]
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

d = Path(sys.argv[1] if len(sys.argv) > 1 else "out_advection")
fig_path = Path(sys.argv[2]) if len(sys.argv) > 2 else d / "advection.png"
plt.style.use("dark_background")


def read_pgm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5"
        w, h = map(int, f.readline().split())
        f.readline()  # 255
        return np.frombuffer(f.read(w * h), dtype=np.uint8).reshape(h, w) / 255.0


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


runs = sorted(p.stem for p in d.glob("*.csv") if p.stem != "summary")
runs.sort(key=lambda r: (r.split("_")[0] != "euler", float(r.split("cfl")[1])))
snaps = [("0", "départ"), ("1", "1/4 tour"), ("2", "1/2 tour"), ("4", "1 tour")]

fig = plt.figure(figsize=(13, 2.2 * len(runs) + 3.2))
gs = fig.add_gridspec(len(runs) + 1, 4, height_ratios=[1] * len(runs) + [1.4])
for r, run in enumerate(runs):
    for c, (k, label) in enumerate(snaps):
        ax = fig.add_subplot(gs[r, c])
        ax.imshow(read_pgm(d / f"{run}_{k}.pgm"), cmap="inferno", vmin=0, vmax=1)
        ax.set_xticks([]), ax.set_yticks([])
        if r == 0:
            ax.set_title(label)
        if c == 0:
            ax.set_ylabel(run.replace("_cfl", "\nCFL "), rotation=0, ha="right", va="center")

colors = {"euler": "#ff7b72", "rk2": "#79c0ff"}
styles = ["-", "--", ":"]
panels = [("max", "max de q"), ("mass_drift", "dérive de masse relative"),
          ("centroid_err_h", "écart du barycentre (en h)"), ("l1", "erreur L1 / solution exacte")]
for c, (key, title) in enumerate(panels):
    ax = fig.add_subplot(gs[-1, c])
    for run in runs:
        scheme, cfl = run.split("_cfl")
        idx = sorted({float(x.split("cfl")[1]) for x in runs}).index(float(cfl))
        data = read_csv(d / f"{run}.csv")
        ax.plot(data["t"], data[key], styles[idx % 3], color=colors[scheme],
                label=f"{scheme} CFL {cfl}")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("t (tours)")
    ax.grid(alpha=0.2)
    if c == 0:
        ax.legend(fontsize=7)

fig.tight_layout()
fig.savefig(fig_path, dpi=110)
print(f"Figure écrite : {fig_path}")
