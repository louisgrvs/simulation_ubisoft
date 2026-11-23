"""Trace les mesures de sand_demo (jour 13) : effondrement d'une colonne, eau contre sable.

Ligne du haut :
    profils finaux h(x)   colonne initiale, eau, sable φ = 20°, 30°, 45°, et sable φ = 30°
                          sans groupes rigides (il s'aplatit : rôle de l'étape 5 de §3.2)
    angle du tas contre φ moyenne sur les graines (barres : écart type), angle des flancs
                          (régression locale) et angle équivalent atan(2H/W) ; droite y = φ
    largeur du dépôt contre φ (même moyenne)
Ligne du bas, au cours du temps (series.csv) :
    largeur du dépôt, énergie cinétique (échelle log), particules par cellule occupée
    (4 au départ : s'il monte, le sable se tasse et perd du volume)

Usage : python tools/plot_sand.py <sorties_sable>
"""
import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

d = Path(sys.argv[1] if len(sys.argv) > 1 else "sorties_sable")


def read_by(path, key):
    """CSV -> {valeur de la colonne key: {colonne: tableau numpy}}, dans l'ordre du fichier."""
    out = OrderedDict()
    if not path.exists():
        return out
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            s = out.setdefault(row[key], {})
            for k, v in row.items():
                if k != key:
                    s.setdefault(k, []).append(float(v))
    return OrderedDict((s, {k: np.array(v) for k, v in cols.items()}) for s, cols in out.items())


plt.style.use("dark_background")
colors = {
    "initial": "#808080",
    "eau": "#7fd1ff",
    "sable 20": "#ffe08a",
    "sable 30": "#ffb86b",
    "sable 45": "#ff6e5a",
    "sable 30 sans groupes": "#c58cff",
    "sable 30 sans 2e proj.": "#8affc1",
}
styles = {"initial": "--", "sable 30 sans groupes": ":", "sable 30 sans 2e proj.": "-."}

fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))
(a1, a2, a3), (a4, a5, a6) = axes

# --- Profils finaux ------------------------------------------------------------------------
profiles = read_by(d / "profiles.csv", "cas")
for name, r in profiles.items():
    c = colors.get(name, "white")
    if name == "initial":
        a1.plot(r["x"], r["h"], color=c, lw=1.5, ls="--", label="colonne initiale")
    elif name == "eau":
        a1.fill_between(r["x"], r["h"], color=c, alpha=0.25)
        a1.plot(r["x"], r["h"], color=c, lw=2, label="eau")
    else:
        a1.plot(r["x"], r["h"], color=c, lw=2, ls=styles.get(name, "-"), label=name.replace("sable", "sable φ ="))
a1.set_title("Dépôt final : hauteur de la plus haute particule")
a1.set_xlabel("x (m)")
a1.set_ylabel("h (m)")
a1.set_xlim(0, 1)
a1.set_ylim(bottom=0)
a1.set_aspect("equal", adjustable="box")

# --- Balayage de φ -------------------------------------------------------------------------
sweep = read_by(d / "sweep.csv", "phi")
if sweep:
    phis = np.array([float(p) for p in sweep])

    def stat(col, keep=None):
        m, s = [], []
        for r in sweep.values():
            v = r[col] if keep is None else r[col][keep(r)]
            m.append(v.mean() if v.size else np.nan)
            s.append(v.std() if v.size else np.nan)
        return np.array(m), np.array(s)

    flank_m, flank_s = stat("angle", keep=lambda r: r["flancs"] > 0)
    eq_m, eq_s = stat("angle_equivalent")
    a2.plot([0, 60], [0, 60], color="gray", lw=1, ls=":", label="angle = φ")
    a2.errorbar(phis, flank_m, yerr=flank_s, color="#ffb86b", marker="o", lw=2, capsize=4,
                label="flancs (régression locale)")
    a2.errorbar(phis, eq_m, yerr=eq_s, color="#7fd1ff", marker="s", mfc="none", lw=1.5, capsize=4,
                label="équivalent atan(2H/W)")
    a2.set_xlim(phis.min() - 3, phis.max() + 3)
    a2.set_ylim(0, 60)
    a2.set_title("Angle du tas contre l'angle de frottement φ")
    a2.set_xlabel("φ (°)")
    a2.set_ylabel("angle mesuré (°)")

    w_m, w_s = stat("largeur")
    a3.errorbar(phis, w_m, yerr=w_s, color="#ffb86b", marker="o", lw=2, capsize=4, label="sable")
    a3.axhline(1.0, color="#7fd1ff", lw=1, ls="--", label="domaine (Lx = 1)")
    a3.set_title("Largeur du dépôt à t_fin")
    a3.set_xlabel("φ (°)")
    a3.set_ylabel("x_max − x_min (m)")
    a3.set_ylim(0, 1.05)
else:
    for ax in (a2, a3):
        ax.text(0.5, 0.5, "sweep.csv absent", ha="center", va="center", transform=ax.transAxes, color="gray")

# --- Séries temporelles --------------------------------------------------------------------
series = read_by(d / "series.csv", "cas")
for name, r in series.items():
    c = colors.get(name, "white")
    ls = styles.get(name, "-")
    lab = name if name == "eau" else name.replace("sable", "sable φ =")
    a4.plot(r["t"], r["largeur"], color=c, lw=2, ls=ls, label=lab)
    a5.semilogy(r["t"][1:], np.maximum(r["Ec"][1:], 1e-12), color=c, lw=1.8, ls=ls, label=lab)
    a6.plot(r["t"], r["particules_par_cellule"], color=c, lw=1.8, ls=ls, label=lab)
a4.set_title("Largeur du dépôt au cours du temps")
a4.set_ylabel("x_max − x_min (m)")
a5.set_title("Énergie cinétique des particules")
a5.set_ylabel("E_c (J/m)")
a6.axhline(4, color="gray", lw=1, ls=":")
a6.set_title("Particules par cellule occupée (4 au départ)")
a6.set_ylabel("particules / cellule")
for ax in (a4, a5, a6):
    ax.set_xlabel("t (s)")

for ax in axes.flat:
    ax.grid(alpha=0.25)
    if ax.lines:
        ax.legend(fontsize=8.5)

fig.suptitle("Jour 13 : effondrement d'une colonne, eau contre sable (Zhu & Bridson §3, PIC)", fontsize=14)
fig.tight_layout()
fig.savefig(d / "sable.png", dpi=130)
print(f"figure : {d / 'sable.png'}")
